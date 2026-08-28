#!/usr/bin/env python3

import csv
import sys
import tempfile
import unittest
from pathlib import Path


SPLICING_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SPLICING_DIR))

import build_neojunction_sequences as spl4  # noqa: E402


class FrameClassificationTests(unittest.TestCase):
    def make_transcript(self, first_length: int, skipped_length: int):
        tx = spl4.Transcript(
            transcript_id="ENSTTEST.1",
            gene_id="ENSGTEST.1",
            gene_name="TEST1",
            gene_type="protein_coding",
            transcript_type="protein_coding",
            chrom="chr1",
            strand="+",
        )
        exon1 = (1, first_length, 1)
        exon2 = (20, 20 + skipped_length - 1, 2)
        exon3 = (40, 51, 3)
        tx.exons_raw = [exon1, exon2, exon3]
        tx.cds_raw = [(start, end, 0) for start, end, _ in tx.exons_raw]
        spl4.finalize_transcript(tx)

        sequences = {
            exon1[:2]: "ATGAAACCC"[:first_length],
            exon2[:2]: "GGGCCCAAAGGG"[:skipped_length],
            exon3[:2]: "CCCGGGAAACCC",
        }
        seq_cache = {
            (tx.chrom, start, end): sequences[(start, end)]
            for start, end, _ in tx.exons_raw
        }
        return tx, seq_cache

    def make_negative_transcript(self, first_length: int, skipped_length: int):
        tx = spl4.Transcript(
            transcript_id="ENSTNEG.1",
            gene_id="ENSGNEG.1",
            gene_name="TESTNEG",
            gene_type="protein_coding",
            transcript_type="protein_coding",
            chrom="chr2",
            strand="-",
        )
        exon1 = (40, 40 + first_length - 1, 1)
        exon2 = (20, 20 + skipped_length - 1, 2)
        exon3 = (1, 12, 3)
        tx.exons_raw = [exon1, exon2, exon3]
        tx.cds_raw = [(start, end, 0) for start, end, _ in tx.exons_raw]
        spl4.finalize_transcript(tx)

        transcript_sequences = {
            exon1[:2]: "ATGAAACCC"[:first_length],
            exon2[:2]: "GGGCCCAAAGGG"[:skipped_length],
            exon3[:2]: "CCCGGGAAACCC",
        }
        seq_cache = {
            (tx.chrom, start, end): spl4.reverse_complement(transcript_sequences[(start, end)])
            for start, end, _ in tx.exons_raw
        }
        return tx, seq_cache

    @staticmethod
    def event_row(tx, acceptor_exon_index: int = 2):
        exon1 = tx.exons_raw[0]
        acceptor_exon = tx.exons_raw[acceptor_exon_index]
        return {
            "sample": "TEST_RNA_TUMOR",
            "chrom": "chr1",
            "strand": "+",
            "star_intron_start": str(exon1[1] + 1),
            "star_intron_end": str(acceptor_exon[0] - 1),
            "junc_id": f"chr1:+:{exon1[1]}-{acceptor_exon[0] - 1}",
            "matched_transcript_id": tx.transcript_id,
            "unique_reads": "20",
            "multimap_reads": "0",
            "total_reads": "25",
            "max_splice_overhang": "60",
            "ssnip_event_class": "exon_skipping",
            "left_feature": "exon_1",
            "right_feature": f"exon_{acceptor_exon_index + 1}",
        }

    def test_inframe_skip_is_not_called_frameshift_when_pipe_splits_codon(self):
        tx, seq_cache = self.make_transcript(first_length=8, skipped_length=9)

        fields = spl4.build_sequences(self.event_row(tx), tx, seq_cache)

        self.assertEqual(fields["cds_length_delta"], "-9")
        self.assertEqual(fields["reading_frame"], "inframe")
        self.assertIn("junction_inside_codon", fields["qc_flags"])

    def test_frameshift_is_detected_when_pipe_lands_on_codon_boundary(self):
        tx, seq_cache = self.make_transcript(first_length=9, skipped_length=8)

        fields = spl4.build_sequences(self.event_row(tx), tx, seq_cache)

        self.assertEqual(fields["cds_length_delta"], "-8")
        self.assertEqual(fields["reading_frame"], "out_of_frame")
        self.assertNotIn("junction_inside_codon", fields["qc_flags"])

    def test_negative_strand_uses_the_same_cds_delta_rule(self):
        tx, seq_cache = self.make_negative_transcript(first_length=8, skipped_length=9)
        first_exon = tx.exons_raw[0]
        last_exon = tx.exons_raw[2]
        row = {
            "chrom": "chr2",
            "strand": "-",
            "star_intron_start": str(last_exon[1] + 1),
            "star_intron_end": str(first_exon[0] - 1),
        }

        fields = spl4.build_sequences(row, tx, seq_cache)

        self.assertEqual(fields["cds_length_delta"], "-9")
        self.assertEqual(fields["reading_frame"], "inframe")
        self.assertTrue(fields["aa_sequence"].startswith("M"))

    def test_unchanged_protein_is_rejected_from_mupexi_output(self):
        tx, seq_cache = self.make_transcript(first_length=9, skipped_length=8)
        row = self.event_row(tx, acceptor_exon_index=1)
        fields = spl4.build_sequences(row, tx, seq_cache)
        confidence, _ = spl4.confidence(row, fields)

        reasons = spl4.sequence_filter_reasons(fields, confidence, False)

        self.assertEqual(fields["reading_frame"], "no_protein_change")
        self.assertIn("frame_no_protein_change", reasons)

    def test_process_file_separates_eligible_and_rejected_events(self):
        tx, seq_cache = self.make_transcript(first_length=8, skipped_length=9)
        eligible = self.event_row(tx)
        unchanged = self.event_row(tx, acceptor_exon_index=1)
        transcripts = {tx.transcript_id: tx, spl4.strip_version(tx.transcript_id): tx}

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_tsv = tmp / "TEST.spl3.event_annotated.tsv"
            fieldnames = list(eligible)
            with input_tsv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows([eligible, unchanged])

            output_tsv = tmp / "TEST.spl4.neojunctions.tsv"
            output_nt = tmp / "TEST.spl4.neojunctions.nt.fa"
            output_aa = tmp / "TEST.spl4.neojunctions.aa.fa"
            rejected_tsv = tmp / "TEST.spl4.sequence_rejected.tsv"
            stats = spl4.process_file(
                input_tsv,
                output_tsv,
                output_nt,
                output_aa,
                rejected_tsv,
                transcripts,
                seq_cache,
                False,
            )

            with output_tsv.open(encoding="utf-8", newline="") as handle:
                output_rows = list(csv.DictReader(handle, delimiter="\t"))
            with rejected_tsv.open(encoding="utf-8", newline="") as handle:
                rejected_rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(stats["output_rows"], 1)
        self.assertEqual(stats["rejected_rows"], 1)
        self.assertEqual(output_rows[0]["reading_frame"], "inframe")
        self.assertEqual(output_rows[0]["sequence_eligible"], "1")
        self.assertIn("frame_no_protein_change", rejected_rows[0]["sequence_filter_reason"])
        self.assertEqual(rejected_rows[0]["peptide_sequence"], "NA")


if __name__ == "__main__":
    unittest.main()
