#!/usr/bin/env python3

import sys
import unittest
from types import SimpleNamespace
from unittest import mock

from splicing import filter_normal_junctions as spl35


def filter_args(max_normal_prevalence):
    return SimpleNamespace(
        max_normal_prevalence=max_normal_prevalence,
        normal_total_samples=0,
        max_normal_sample_count=None,
    )


class NormalFilteringTests(unittest.TestCase):
    def test_cli_default_is_strict(self):
        with mock.patch.object(
            sys,
            "argv",
            ["filter_normal_junctions.py", "--root", ".", "--normal-ref", "normal.tsv"],
        ):
            args = spl35.parse_args()
        self.assertEqual(args.max_normal_prevalence, 0.0)

    def test_strict_mode_keeps_only_normal_absent_rows(self):
        absent, keep_absent = spl35.annotate_row({}, None, filter_args(0.0))
        hit = spl35.NormalHit()
        hit.sample_count = 1
        hit.prevalence = 1 / 9166
        present, keep_present = spl35.annotate_row({}, hit, filter_args(0.0))

        self.assertTrue(keep_absent)
        self.assertEqual(absent["normal_status"], "absent")
        self.assertFalse(keep_present)
        self.assertEqual(present["normal_status"], "filtered_normal_present")

    def test_permissive_mode_retains_rare_normal_match(self):
        hit = spl35.NormalHit()
        hit.sample_count = 7
        hit.prevalence = 0.0007637

        annotated, keep = spl35.annotate_row({}, hit, filter_args(0.01))

        self.assertTrue(keep)
        self.assertEqual(annotated["normal_status"], "present_below_threshold")


if __name__ == "__main__":
    unittest.main()
