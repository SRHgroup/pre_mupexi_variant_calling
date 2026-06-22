#!/usr/bin/env python3
"""Lift over compact normal splice-junction references between genome builds."""

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
import tempfile
from bisect import bisect_left
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple


EXTRA_FIELDS = [
    "liftover_status",
    "liftover_from_build",
    "liftover_to_build",
    "source_chrom",
    "source_left_boundary",
    "source_right_boundary",
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Lift a compact normal splice-junction reference by independently "
            "lifting the left and right 1-bp splice boundary positions."
        )
    )
    ap.add_argument("--input", required=True, help="Input compact normal reference TSV/TSV.GZ")
    ap.add_argument("--out", required=True, help="Lifted output TSV/TSV.GZ")
    ap.add_argument("--chain", required=True, help="UCSC chain file, e.g. hg19ToHg38.over.chain.gz")
    ap.add_argument(
        "--engine",
        choices=("python", "ucsc"),
        default="python",
        help="Liftover engine. python parses the chain file directly; ucsc shells out to UCSC liftOver.",
    )
    ap.add_argument("--liftover-bin", default="liftOver", help="UCSC liftOver executable, only used with --engine ucsc")
    ap.add_argument("--from-build", default="GRCh37", help="Label written to liftover_from_build")
    ap.add_argument("--to-build", default="GRCh38", help="Label written to liftover_to_build")
    ap.add_argument(
        "--unmapped-out",
        default="",
        help="Optional unmapped audit TSV. Default: OUT with .unmapped.tsv appended before .gz.",
    )
    ap.add_argument("--keep-cross-chrom", action="store_true", help="Keep junctions whose two boundaries lift to different chromosomes")
    ap.add_argument(
        "--allow-inverted",
        action="store_true",
        help="Keep junctions where lifted left_boundary is greater than lifted right_boundary",
    )
    ap.add_argument("--tmpdir", default="", help="Optional temporary directory")
    return ap.parse_args()


def open_input(path: Path):
    if path.suffix in {".gz", ".bgz"}:
        return gzip.open(path, "rt")
    return path.open("r", encoding="utf-8", newline="")


def open_output(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        return gzip.open(path, "wt")
    return path.open("w", encoding="utf-8", newline="")


def default_unmapped_path(out_path: Path, explicit: str) -> Path:
    if explicit:
        return Path(explicit)
    if out_path.name.endswith(".gz"):
        return out_path.with_name(f"{out_path.name[:-3]}.unmapped.tsv")
    return out_path.with_name(f"{out_path.name}.unmapped.tsv")


def parse_int(value: str, label: str) -> int:
    try:
        return int(str(value).strip())
    except ValueError as exc:
        raise SystemExit(f"ERROR: cannot parse integer {label}: {value!r}") from exc


def chrom_value(row: Dict[str, str]) -> str:
    chrom = row.get("chrom", "")
    if not chrom:
        raise SystemExit("ERROR: input reference is missing required column: chrom")
    return chrom


def boundary_values(row: Dict[str, str]) -> Tuple[int, int]:
    left = row.get("left_boundary", "")
    right = row.get("right_boundary", "")
    if not left or not right:
        raise SystemExit("ERROR: input reference must contain left_boundary and right_boundary columns")
    return parse_int(left, "left_boundary"), parse_int(right, "right_boundary")


def point_to_bed(chrom: str, pos_1based: int, name: str) -> str:
    if pos_1based < 1:
        raise SystemExit(f"ERROR: cannot liftover non-positive 1-based position {pos_1based} for {name}")
    return f"{chrom}\t{pos_1based - 1}\t{pos_1based}\t{name}\n"


def open_chain(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r", encoding="utf-8", newline="")


def run_liftover(liftover_bin: str, bed_in: Path, chain: Path, bed_out: Path, unmapped_bed: Path) -> None:
    binary = shutil.which(liftover_bin)
    if binary is None:
        binary = shutil.which("liftOver") if liftover_bin != "liftOver" else None
    if binary is None:
        raise SystemExit(
            f"ERROR: liftOver executable not found: {liftover_bin}. "
            "Load a UCSC/liftOver module or pass --liftover-bin /path/to/liftOver."
        )
    cmd = [binary, str(bed_in), str(chain), str(bed_out), str(unmapped_bed)]
    print(f"[liftover-normal-ref] running: {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True)


def split_boundary_name(name: str) -> Tuple[int, str]:
    try:
        row_index, side = name.split("|", 1)
        if side not in {"L", "R"}:
            raise ValueError
        return int(row_index), side
    except ValueError as exc:
        raise SystemExit(f"ERROR: malformed temporary BED name from liftOver: {name!r}") from exc


def read_lifted_bed(path: Path) -> Dict[int, Dict[str, Tuple[str, int, int]]]:
    lifted: Dict[int, Dict[str, Tuple[str, int, int]]] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start0, end0, name = parts[:4]
            row_index, side = split_boundary_name(name)
            lifted.setdefault(row_index, {})[side] = (chrom, int(start0), int(end0))
    return lifted


def read_unmapped_names(path: Path) -> Dict[int, set]:
    names: Dict[int, set] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                row_index, side = split_boundary_name(parts[3])
                names.setdefault(row_index, set()).add(side)
    return names


def collect_boundary_positions(input_path: Path) -> Tuple[List[str], int, Dict[str, List[Tuple[int, int, str]]]]:
    positions: Dict[str, List[Tuple[int, int, str]]] = {}
    with open_input(input_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"ERROR: input reference has no header: {input_path}")
        fieldnames = reader.fieldnames
        input_rows = 0
        for idx, row in enumerate(reader):
            chrom = chrom_value(row)
            left, right = boundary_values(row)
            if left < 1 or right < 1:
                raise SystemExit(f"ERROR: cannot liftover non-positive boundary in row {idx + 2}")
            positions.setdefault(chrom, []).append((left - 1, idx, "L"))
            positions.setdefault(chrom, []).append((right - 1, idx, "R"))
            input_rows += 1
    return fieldnames, input_rows, positions


def map_positions_in_block(
    entries: List[Tuple[int, int, str]],
    coords: List[int],
    lifted: Dict[int, Dict[str, Tuple[str, int, int]]],
    source_start0: int,
    source_end0: int,
    target_chrom: str,
    target_start0: int,
    source_strand: str,
) -> int:
    mapped = 0
    lo = bisect_left(coords, source_start0)
    hi = bisect_left(coords, source_end0)
    for source_pos0, row_index, side in entries[lo:hi]:
        lifted_row = lifted.setdefault(row_index, {})
        if side in lifted_row:
            continue
        if source_strand == "+":
            target_pos0 = target_start0 + (source_pos0 - source_start0)
        else:
            target_pos0 = target_start0 + ((source_end0 - 1) - source_pos0)
        lifted_row[side] = (target_chrom, target_pos0, target_pos0 + 1)
        mapped += 1
    return mapped


def lift_positions_with_python(
    input_path: Path, chain_path: Path
) -> Tuple[List[str], int, Dict[int, Dict[str, Tuple[str, int, int]]], Dict[int, Set[str]], int]:
    fieldnames, input_rows, positions = collect_boundary_positions(input_path)
    for chrom in positions:
        positions[chrom].sort(key=lambda item: item[0])
    coords_by_chrom = {chrom: [entry[0] for entry in entries] for chrom, entries in positions.items()}
    lifted: Dict[int, Dict[str, Tuple[str, int, int]]] = {}
    mapped_boundaries = 0

    active = None
    t_cursor = 0
    q_cursor = 0
    with open_chain(chain_path) as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                active = None
                continue
            parts = line.split()
            if parts[0] == "chain":
                if len(parts) < 13:
                    raise SystemExit(f"ERROR: malformed chain header: {line}")
                t_chrom = parts[2]
                t_strand = parts[4]
                t_start = parse_int(parts[5], "chain tStart")
                q_chrom = parts[7]
                q_size = parse_int(parts[8], "chain qSize")
                q_strand = parts[9]
                q_start = parse_int(parts[10], "chain qStart")
                if t_strand != "+":
                    active = None
                    continue
                if q_strand not in {"+", "-"}:
                    active = None
                    continue
                active = (t_chrom, q_chrom, q_size, q_strand)
                t_cursor = t_start
                q_cursor = q_start
                continue
            if active is None:
                continue

            block = [parse_int(value, "chain block") for value in parts]
            if len(block) not in {1, 3}:
                raise SystemExit(f"ERROR: malformed chain block: {line}")
            size = block[0]
            t_chrom, q_chrom, q_size, q_strand = active
            entries = positions.get(q_chrom)
            coords = coords_by_chrom.get(q_chrom)
            if entries and coords:
                if q_strand == "+":
                    source_start0 = q_cursor
                    source_end0 = q_cursor + size
                else:
                    source_start0 = q_size - (q_cursor + size)
                    source_end0 = q_size - q_cursor
                mapped_boundaries += map_positions_in_block(
                    entries,
                    coords,
                    lifted,
                    source_start0,
                    source_end0,
                    t_chrom,
                    t_cursor,
                    q_strand,
                )
            t_cursor += size
            q_cursor += size
            if len(block) == 3:
                t_cursor += block[1]
                q_cursor += block[2]

    return fieldnames, input_rows, lifted, {}, mapped_boundaries


def lift_positions_with_ucsc(
    input_path: Path, chain_path: Path, liftover_bin: str, tmpdir_arg: str
) -> Tuple[List[str], int, Dict[int, Dict[str, Tuple[str, int, int]]], Dict[int, set], int]:
    with open_input(input_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"ERROR: input reference has no header: {input_path}")
        fieldnames = reader.fieldnames
        input_rows = 0
        with tempfile.TemporaryDirectory(dir=tmpdir_arg or None) as tmp:
            tmpdir = Path(tmp)
            bed_in = tmpdir / "normal_reference.boundaries.bed"
            bed_out = tmpdir / "normal_reference.boundaries.lifted.bed"
            bed_unmapped = tmpdir / "normal_reference.boundaries.unmapped.bed"
            with bed_in.open("w", encoding="utf-8", newline="") as bed:
                for idx, row in enumerate(reader):
                    chrom = chrom_value(row)
                    left, right = boundary_values(row)
                    bed.write(point_to_bed(chrom, left, f"{idx}|L"))
                    bed.write(point_to_bed(chrom, right, f"{idx}|R"))
                    input_rows += 1

            run_liftover(liftover_bin, bed_in, chain_path, bed_out, bed_unmapped)
            lifted = read_lifted_bed(bed_out)
            unmapped_names = read_unmapped_names(bed_unmapped)
    mapped_boundaries = sum(len(value) for value in lifted.values())
    return fieldnames, input_rows, lifted, unmapped_names, mapped_boundaries


def append_fields(fieldnames: Sequence[str], extra: Sequence[str]) -> List[str]:
    fields = list(fieldnames)
    for field in extra:
        if field not in fields:
            fields.append(field)
    return fields


def lift_reference(args: argparse.Namespace) -> Dict[str, int]:
    input_path = Path(args.input)
    out_path = Path(args.out)
    chain_path = Path(args.chain)
    unmapped_path = default_unmapped_path(out_path, args.unmapped_out)
    if not input_path.exists():
        raise SystemExit(f"ERROR: input reference does not exist: {input_path}")
    if not chain_path.exists():
        raise SystemExit(f"ERROR: chain file does not exist: {chain_path}")

    if args.engine == "ucsc":
        fieldnames, input_rows, lifted, unmapped_names, mapped_boundaries = lift_positions_with_ucsc(
            input_path, chain_path, args.liftover_bin, args.tmpdir
        )
    else:
        print(f"[liftover-normal-ref] using pure-Python chain parser: {chain_path}", file=sys.stderr)
        fieldnames, input_rows, lifted, unmapped_names, mapped_boundaries = lift_positions_with_python(input_path, chain_path)

    stats = {
        "input_rows": input_rows,
        "output_rows": 0,
        "unmapped_rows": 0,
        "cross_chrom_rows": 0,
        "inverted_rows": 0,
        "partially_unmapped_rows": 0,
        "mapped_boundaries": mapped_boundaries,
    }
    out_fields = append_fields(fieldnames, EXTRA_FIELDS)
    unmapped_fields = append_fields(fieldnames, EXTRA_FIELDS + ["liftover_reason"])

    unmapped_path.parent.mkdir(parents=True, exist_ok=True)
    with open_output(out_path) as out, unmapped_path.open("w", encoding="utf-8", newline="") as unmapped:
        writer = csv.DictWriter(out, fieldnames=out_fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        unmapped_writer = csv.DictWriter(
            unmapped, fieldnames=unmapped_fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        unmapped_writer.writeheader()

        with open_input(input_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for idx, row in enumerate(reader):
                source_chrom = chrom_value(row)
                source_left, source_right = boundary_values(row)
                lifted_row = lifted.get(idx, {})
                left = lifted_row.get("L")
                right = lifted_row.get("R")
                output_row = dict(row)
                output_row.update(
                    {
                        "liftover_from_build": args.from_build,
                        "liftover_to_build": args.to_build,
                        "source_chrom": source_chrom,
                        "source_left_boundary": str(source_left),
                        "source_right_boundary": str(source_right),
                    }
                )

                if left is None or right is None:
                    stats["unmapped_rows"] += 1
                    if left is None and right is None:
                        reason = "both_boundaries_unmapped"
                    else:
                        stats["partially_unmapped_rows"] += 1
                        reason = "one_boundary_unmapped"
                    if args.engine == "ucsc" and idx not in unmapped_names and not lifted_row:
                        reason = "both_boundaries_missing_from_liftover_output"
                    output_row.update({"liftover_status": "unmapped", "liftover_reason": reason})
                    unmapped_writer.writerow(output_row)
                    continue

                left_chrom, left_start0, left_end0 = left
                right_chrom, right_start0, right_end0 = right
                left_pos = left_start0 + 1
                right_pos = right_start0 + 1
                if left_end0 - left_start0 != 1 or right_end0 - right_start0 != 1:
                    stats["unmapped_rows"] += 1
                    output_row.update({"liftover_status": "unmapped", "liftover_reason": "boundary_not_1bp_after_liftover"})
                    unmapped_writer.writerow(output_row)
                    continue
                if left_chrom != right_chrom and not args.keep_cross_chrom:
                    stats["cross_chrom_rows"] += 1
                    stats["unmapped_rows"] += 1
                    output_row.update({"liftover_status": "unmapped", "liftover_reason": "boundaries_lift_to_different_chroms"})
                    unmapped_writer.writerow(output_row)
                    continue
                if left_pos > right_pos and not args.allow_inverted:
                    stats["inverted_rows"] += 1
                    stats["unmapped_rows"] += 1
                    output_row.update({"liftover_status": "unmapped", "liftover_reason": "lifted_left_greater_than_right"})
                    unmapped_writer.writerow(output_row)
                    continue

                output_row.update(
                    {
                        "chrom": left_chrom,
                        "left_boundary": str(left_pos),
                        "right_boundary": str(right_pos),
                        "liftover_status": "lifted",
                    }
                )
                writer.writerow(output_row)
                stats["output_rows"] += 1

    summary_path = out_path.with_name(f"{out_path.name[:-3]}.summary.tsv" if out_path.name.endswith(".gz") else f"{out_path.name}.summary.tsv")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8", newline="") as summary:
        summary_writer = csv.writer(summary, delimiter="\t", lineterminator="\n")
        summary_writer.writerow(["metric", "value"])
        summary_writer.writerow(["input", str(input_path)])
        summary_writer.writerow(["output", str(out_path)])
        summary_writer.writerow(["unmapped_output", str(unmapped_path)])
        summary_writer.writerow(["chain", str(chain_path)])
        summary_writer.writerow(["engine", args.engine])
        for key, value in stats.items():
            summary_writer.writerow([key, value])

    print(
        "[liftover-normal-ref][done] "
        + " ".join(f"{key}={value}" for key, value in stats.items()),
        file=sys.stderr,
    )
    print(f"[liftover-normal-ref] unmapped: {unmapped_path}", file=sys.stderr)
    print(f"[liftover-normal-ref] summary: {summary_path}", file=sys.stderr)
    return stats


def main() -> int:
    csv.field_size_limit(sys.maxsize)
    lift_reference(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
