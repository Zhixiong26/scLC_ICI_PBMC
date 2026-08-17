from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from mvi_utils import (  # noqa: E402
    aggregate_allc,
    canonical_cell_id,
    infer_sample_id,
    load_annotations,
    load_sample_metadata,
    region_lookup,
    regions_from_var,
)


class PipelineUtilsTests(unittest.TestCase):
    def test_canonical_cell_id_only_changes_first_delimiter(self):
        self.assertEqual(canonical_cell_id("D01-AA-BB"), "D01_AA-BB")
        self.assertEqual(canonical_cell_id("D01_AA-BB"), "D01_AA-BB")

    def test_scanpy_and_allcools_cell_ids_match(self):
        expected = "IR01__AAAGAAGAAGAGGAGAG"
        self.assertEqual(canonical_cell_id("IR01_AAAGAAGAAGAGGAGAG"), expected)
        self.assertEqual(canonical_cell_id("IR01__AAAGAAGAAGAGGAGAG"), expected)
        self.assertEqual(
            canonical_cell_id("25110891_IR01_Met__AAAGAAGAAGAGGAGAG"),
            expected,
        )

    def test_methscan_sample_id_regex(self):
        pattern = r"^(?:[^_]+_)?((?:IR|NR)[0-9][0-9])(?:_Met)?(?:__|_)"
        self.assertEqual(infer_sample_id("25110891_IR01_Met__barcode", pattern), "IR01")
        self.assertEqual(infer_sample_id("NR05__barcode", pattern), "NR05")

    def test_parse_allcools_region_names(self):
        var = pd.DataFrame(index=["chr1-0-5000", "chr1-10000-15000"])
        regions = regions_from_var(var)
        self.assertEqual(regions["bin"].tolist(), [0, 2])

    def test_parse_real_allcools_coordinate_columns(self):
        var = pd.DataFrame(
            {"chrom": ["chr1", "chr1"], "start": [0, 5000], "end": [5000, 10000]},
            index=["chr1_0", "chr1_1"],
        )
        regions = regions_from_var(var)
        self.assertEqual(regions["feature_index"].tolist(), [0, 1])

    def test_aggregate_selected_cgn_sites(self):
        var = pd.DataFrame(index=["chr1-0-5000", "chr1-5000-10000"])
        lookup = region_lookup(regions_from_var(var))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cell.allc.tsv.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("chr1\t1\t+\tCGN\t1\t2\t1\n")
                handle.write("chr1\t5000\t+\tCGN\t2\t3\t1\n")
                handle.write("chr1\t5001\t+\tCGN\t4\t5\t1\n")
                handle.write("chr1\t6000\t+\tCHN\t8\t9\t1\n")
            indices, mc, cov, stats = aggregate_allc(path, lookup, 2)
        self.assertEqual(indices.tolist(), [0, 1])
        self.assertEqual(mc.tolist(), [3, 4])
        self.assertEqual(cov.tolist(), [5, 5])
        self.assertEqual(stats["selected_sites"], 3)

    def test_sample_metadata_and_annotation_merge(self):
        with tempfile.TemporaryDirectory() as directory:
            metadata_path = Path(directory) / "sample_metadata.tsv"
            metadata_path.write_text("sample_id\tcondition\nIR_01\tIR\nNR_01\tNR\n")
            metadata = load_sample_metadata(metadata_path)
            self.assertEqual(metadata.loc["IR_01", "condition"], "IR")
            annotations, stats = load_annotations(
                ["IR_01_cellA", "NR_01_cellB"], None, metadata_path
            )
            self.assertEqual(annotations["sample_id"].tolist(), ["IR_01", "NR_01"])
            self.assertEqual(annotations["condition"].tolist(), ["IR", "NR"])
            self.assertEqual(stats["sample_metadata_rows"], 2)

    def test_scanpy_annotation_column_aliases(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata_path = root / "sample_metadata.tsv"
            metadata_path.write_text("sample_id\tcondition\nIR01\tIR\n")
            annotation_path = root / "scanpy_annotation.csv"
            annotation_path.write_text(
                "cell_id,sample,group,cell_type_integrated\n"
                "IR01_AAAC,IR01,IR,T_cell\n"
            )
            annotations, stats = load_annotations(
                ["IR01__AAAC"],
                annotation_path,
                metadata_path,
                r"^((?:IR|NR)[0-9]{2})(?:__|_)",
            )
            self.assertEqual(annotations.loc["IR01__AAAC", "sample_id"], "IR01")
            self.assertEqual(annotations.loc["IR01__AAAC", "condition"], "IR")
            self.assertEqual(annotations.loc["IR01__AAAC", "cell_type"], "T_cell")
            self.assertEqual(stats["fully_annotated_selected_cells"], 1)
            self.assertEqual(stats["cell_type_annotated_selected_cells"], 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
