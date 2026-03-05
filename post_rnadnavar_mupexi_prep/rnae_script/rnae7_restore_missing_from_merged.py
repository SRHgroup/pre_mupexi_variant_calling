#!/usr/bin/env python3
import argparse
import gzip
import sys


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_key(cols):
    return (cols[0], int(cols[1]), cols[3], cols[4])


def cmp_key(a, b):
    if a[0] != b[0]:
        return -1 if a[0] < b[0] else 1
    if a[1] != b[1]:
        return -1 if a[1] < b[1] else 1
    if a[2] != b[2]:
        return -1 if a[2] < b[2] else 1
    if a[3] != b[3]:
        return -1 if a[3] < b[3] else 1
    return 0


def gt_missing(sample_value, fmt):
    keys = fmt.split(":")
    try:
        gt_idx = keys.index("GT")
    except ValueError:
        return True
    vals = sample_value.split(":")
    gt = vals[gt_idx] if gt_idx < len(vals) else "."
    return gt in (".", "./.", ".|.", "")


def read_header_and_samples(path):
    fh = open_text(path, "rt")
    header = []
    samples = []
    for line in fh:
        header.append(line)
        if line.startswith("#CHROM"):
            cols = line.rstrip("\n").split("\t")
            samples = cols[9:]
            break
    return fh, header, samples


def next_variant_line(fh):
    for line in fh:
        if not line.startswith("#"):
            return line
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--merged", required=True)
    ap.add_argument("--genotyped", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--normal-name", required=True)
    ap.add_argument("--dna-name", required=True)
    ap.add_argument("--rna-name", required=True)
    args = ap.parse_args()

    mfh, _mheader, msamples = read_header_and_samples(args.merged)
    gfh, gheader, gsamples = read_header_and_samples(args.genotyped)

    want = [args.normal_name, args.dna_name, args.rna_name]
    m_idx = {}
    g_idx = {}
    for s in want:
        if s in msamples:
            m_idx[s] = 9 + msamples.index(s)
        if s in gsamples:
            g_idx[s] = 9 + gsamples.index(s)

    with open_text(args.out, "wt") as out:
        for h in gheader:
            out.write(h)

        mline = next_variant_line(mfh)
        restored = 0
        total = 0
        while True:
            gline = next_variant_line(gfh)
            if gline is None:
                break
            total += 1
            gcols = gline.rstrip("\n").split("\t")
            gkey = parse_key(gcols)

            mcols = None
            while mline is not None:
                cand = mline.rstrip("\n").split("\t")
                mkey = parse_key(cand)
                c = cmp_key(mkey, gkey)
                if c < 0:
                    mline = next_variant_line(mfh)
                    continue
                if c == 0:
                    mcols = cand
                break

            if mcols is not None:
                gfmt = gcols[8]
                mfmt = mcols[8]
                for s in want:
                    if s not in g_idx or s not in m_idx:
                        continue
                    gi = g_idx[s]
                    mi = m_idx[s]
                    if gi >= len(gcols) or mi >= len(mcols):
                        continue
                    if gt_missing(gcols[gi], gfmt) and not gt_missing(mcols[mi], mfmt):
                        gcols[gi] = mcols[mi]
                        restored += 1

            out.write("\t".join(gcols) + "\n")

    print(f"[info] restore-missing: restored_sample_fields={restored} total_records={total}", file=sys.stderr)


if __name__ == "__main__":
    main()

