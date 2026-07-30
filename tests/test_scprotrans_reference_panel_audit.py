#!/usr/bin/env python3
"""Tests for the scProTrans reference-panel feasibility audit."""

import gzip
import importlib.util
import os
import subprocess
import tempfile
import unittest

import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT_PATH = os.path.normpath(
    os.path.join(HERE, "..", "scripts", "38_scprotrans_reference_panel_audit.py")
)


def load_module():
    spec = importlib.util.spec_from_file_location("scprotrans_reference_audit", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestProteinSemantics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_module()

    def test_hla_dr_labels_map_to_complex_not_single_chain(self):
        for label in ("HLA.DR", "HLA-DR-A0159", "MHCII(HLA-DR)"):
            self.assertEqual(self.mod.canonical_target(label), "HLA_DR_COMPLEX")

        spec = self.mod.TARGET_SPECS["HLA_DR_COMPLEX"]
        self.assertEqual(spec["adt_resolution"], "heterodimer_or_complex_adt")
        self.assertEqual(spec["sequence_strategy"], "mean_chain_embedding")
        self.assertEqual(spec["query_genes"], ("HLA-DRA", "HLA-DRB1"))

    def test_pan_hla_ii_is_not_collapsed_to_hla_dr(self):
        self.assertEqual(
            self.mod.canonical_target("HLA-DR-DP-DQ-A1018"),
            "HLA_II_PAN",
        )

    def test_priority_targets_are_normalized(self):
        expected = {
            "CD74-A0935": "CD74",
            "PD-L1(CD274)": "PD_L1",
            "CD274-A0007": "PD_L1",
            "EpCAM(CD326)": "EPCAM",
            "CD47-A0026": "CD47",
            "HLA-A-B-C-A0058": "HLA_I_PAN",
            "B2M": "B2M",
        }
        observed = {label: self.mod.canonical_target(label) for label in expected}
        self.assertEqual(observed, expected)


class TestStreamingPanelReaders(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_module()

    def test_csv_row_reader_skips_blank_index_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "adt.csv.gz")
            with gzip.open(path, "wt", newline="") as handle:
                handle.write(",cell1,cell2\n")
                handle.write("HLA-DR-A0159,1,2\n")
                handle.write("CD74-A0935,3,4\n")

            labels, cell_count = self.mod.read_row_feature_labels(path, delimiter=",")

        self.assertEqual(labels, ["HLA-DR-A0159", "CD74-A0935"])
        self.assertEqual(cell_count, 2)

    def test_tsv_row_reader_skips_nonblank_cell_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "adt.tsv.gz")
            with gzip.open(path, "wt", newline="") as handle:
                handle.write("cell1\tcell2\n")
                handle.write("HLA.DR\t1\t2\n")
                handle.write("CD47\t3\t4\n")

            labels, cell_count = self.mod.read_row_feature_labels(path, delimiter="\t")

        self.assertEqual(labels, ["HLA.DR", "CD47"])
        self.assertEqual(cell_count, 2)

    def test_panel_summary_counts_measured_targets(self):
        labels = [
            "HLA-DR-A0159",
            "HLA-DR-DP-DQ-A1018",
            "CD74-A0935",
            "CD274-A0007",
            "Mouse-IgG1-k-Ctrl-A0090",
        ]
        summary = self.mod.summarize_panel(labels)
        self.assertEqual(summary["adt_count"], 5)
        self.assertEqual(summary["target_adt_count"], 4)
        self.assertTrue(summary["has_hla_dr_mhcii"])
        self.assertTrue(summary["has_cd74"])
        self.assertTrue(summary["has_pd_l1"])

    def test_row_matrix_index_returns_feature_and_cell_names(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "adt.tsv.gz")
            with gzip.open(path, "wt", newline="") as handle:
                handle.write("\tbarcode\tcell_b\tcell_a\n")
                handle.write("HLA-DR.ADT\t1\t2\t3\n")
                handle.write("CD274_PD-L1\t4\t5\t6\n")

            labels, cells = self.mod.read_row_feature_matrix_index(path, delimiter="\t")

        self.assertEqual(labels, ["HLA-DR.ADT", "CD274_PD-L1"])
        self.assertEqual(cells, ["barcode", "cell_b", "cell_a"])


class TestStructuredReferenceReaders(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_module()

    def test_old_h5ad_categorical_codes_are_decoded(self):
        import h5py

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "old_encoding.h5ad")
            with h5py.File(path, "w") as handle:
                obs = handle.create_group("obs")
                obs.create_dataset("_index", data=np.asarray([b"c1", b"c2", b"c3"]))
                obs.create_dataset("DonorNumber", data=np.asarray([0, 1, 0], dtype=np.int8))
                obs.create_dataset("batch", data=np.asarray([0, 1, 1], dtype=np.int8))
                obs_categories = obs.create_group("__categories")
                obs_categories.create_dataset("DonorNumber", data=np.asarray([b"d1", b"d2"]))
                obs_categories.create_dataset("batch", data=np.asarray([b"b1", b"b2"]))

                var = handle.create_group("var")
                var.create_dataset(
                    "_index", data=np.asarray([b"GAPDH", b"HLA-DR", b"CD274"])
                )
                var.create_dataset("feature_types", data=np.asarray([0, 1, 1], dtype=np.int8))
                var_categories = var.create_group("__categories")
                var_categories.create_dataset(
                    "feature_types", data=np.asarray([b"GEX", b"ADT"])
                )

            labels, cells, donors, batches = self.mod.read_h5ad_panel(path)

        self.assertEqual(labels, ["HLA-DR", "CD274"])
        self.assertEqual(cells, 3)
        self.assertEqual(donors, 2)
        self.assertEqual(batches, 2)

    def test_appledouble_files_are_excluded(self):
        paths = [
            "/tmp/._GSM8061744_ADT.txt.gz",
            "/tmp/GSM8061744_ADT.txt.gz",
            "/tmp/.DS_Store",
        ]
        self.assertEqual(
            self.mod.exclude_metadata_sidecars(paths),
            ["/tmp/GSM8061744_ADT.txt.gz"],
        )

    def test_double_gzip_rds_reader_recovers_dimensions_and_barcodes(self):
        if not self.mod.rscript_available():
            self.skipTest("Rscript is required for the double-gzip RDS contract")

        with tempfile.TemporaryDirectory() as tmpdir:
            inner = os.path.join(tmpdir, "matrix.rds.gz")
            outer = os.path.join(tmpdir, "matrix.rds.gz.gz")
            r_code = (
                "x <- matrix(1:6, nrow=2, dimnames=list(c('G1','G2'), "
                "c('cell_a','cell_b','cell_c'))); "
                "con <- gzfile(commandArgs(TRUE)[1], 'wb'); saveRDS(x, con); close(con)"
            )
            subprocess.run(
                ["Rscript", "-e", r_code, inner],
                check=True,
                capture_output=True,
                text=True,
            )
            with open(outer, "wb") as output:
                subprocess.run(["gzip", "-c", inner], check=True, stdout=output)

            rows, columns, cells = self.mod.read_double_gzip_rds_index(outer)

        self.assertEqual((rows, columns), (2, 3))
        self.assertEqual(cells, {"cell_a", "cell_b", "cell_c"})

    def test_multimodal_alignment_is_intersection_based(self):
        summary = self.mod.summarize_cell_alignment(
            rna_cells={"a", "b", "c", "rna_only"},
            adt_cells={"a", "b", "c", "adt_only"},
            csp_cells={"a", "b", "c", "csp_only"},
            hto_cells={"a", "b", "c", "hto_only"},
        )

        self.assertEqual(summary["complete_cell_count"], 3)
        self.assertEqual(summary["rna_only_count"], 1)
        self.assertEqual(summary["adt_only_count"], 1)
        self.assertEqual(summary["csp_only_count"], 1)
        self.assertEqual(summary["hto_only_count"], 1)

    def test_epithelial_bridge_arm_audits_real_multimodal_overlap(self):
        if not self.mod.rscript_available():
            self.skipTest("Rscript is required for the epithelial bridge audit")

        with tempfile.TemporaryDirectory() as tmpdir:
            rna_path = os.path.join(tmpdir, "rna.rds.gz")
            inner_path = os.path.join(tmpdir, "inner.rds.gz")
            r_code = (
                "x <- matrix(1:8, nrow=2, dimnames=list(c('HLA-DRA','CD274'), "
                "c('a','b','c','rna_only'))); "
                "con <- gzfile(commandArgs(TRUE)[1], 'wb'); saveRDS(x, con); close(con)"
            )
            subprocess.run(
                ["Rscript", "-e", r_code, inner_path],
                check=True,
                capture_output=True,
                text=True,
            )
            with open(rna_path, "wb") as output:
                subprocess.run(["gzip", "-c", inner_path], check=True, stdout=output)

            modality_specs = {
                "adt.tsv.gz": (["HLA-DR.ADT"], ["a", "b", "c", "adt_only"]),
                "csp.tsv.gz": (["CD274_PD-L1"], ["a", "b", "c", "csp_only"]),
                "hto.tsv.gz": (["Hashtag 1"], ["a", "b", "c", "hto_only"]),
            }
            for filename, (features, cells) in modality_specs.items():
                with gzip.open(os.path.join(tmpdir, filename), "wt", newline="") as handle:
                    handle.write("\t" + "\t".join(cells) + "\n")
                    for feature in features:
                        handle.write(feature + "\t" + "\t".join("1" for _ in cells) + "\n")

            row, targets = self.mod.audit_gse254985_arm(
                dataset_arm="GSE254985_HP20276_untreated",
                accession="GSE254985/GSM-test",
                treatment="untreated",
                rna_path=rna_path,
                adt_path=os.path.join(tmpdir, "adt.tsv.gz"),
                csp_path=os.path.join(tmpdir, "csp.tsv.gz"),
                hto_path=os.path.join(tmpdir, "hto.tsv.gz"),
            )

        self.assertEqual(row["biological_context"], "human_pancreatic_islet_epithelial")
        self.assertEqual(row["cell_count"], 3)
        self.assertEqual(row["rna_feature_count"], 2)
        self.assertEqual(row["rna_cell_count"], 4)
        self.assertEqual(row["complete_cell_count"], 3)
        self.assertEqual(targets, {"HLA_DR_COMPLEX", "PD_L1"})

    def test_tsv_writer_unions_fields_across_dataset_types(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "audit.tsv")
            self.mod._write_tsv(
                path=self.mod.Path(path),
                rows=[{"dataset": "immune"}, {"dataset": "epithelial", "rna_cells": 3}],
            )
            with open(path, "r", encoding="utf-8") as handle:
                lines = handle.read().splitlines()

        self.assertEqual(lines[0], "dataset\trna_cells")
        self.assertEqual(lines[1], "immune\t")
        self.assertEqual(lines[2], "epithelial\t3")


class TestEmbeddingAndEligibilityAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_module()

    def test_npz_gene_embedding_reader_reports_dimension(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "gene_embedding.npz")
            np.savez(
                path,
                gene=np.array(["HLA-DRA", "HLA-DRB1", "CD74"]),
                embedding=np.zeros((3, 100), dtype=np.float32),
            )
            genes, dimension = self.mod.read_gene_embedding_index(path)

        self.assertEqual(genes, {"HLA-DRA", "HLA-DRB1", "CD74"})
        self.assertEqual(dimension, 100)

    def test_measured_hla_dr_with_complete_embeddings_is_benchmark_eligible(self):
        measured = {
            "GSE164378_3P": {"HLA_DR_COMPLEX", "PD_L1", "CD47"},
            "GSE128639_BMMC": {"HLA_DR_COMPLEX"},
        }
        rows = self.mod.build_target_coverage_rows(
            measured_targets=measured,
            protein_entries={"P01903", "P01911", "Q9NZQ7", "Q08722"},
            gene_names={"HLA-DRA", "HLA-DRB1", "CD274", "CD47"},
        )
        hla = next(row for row in rows if row["target_id"] == "HLA_DR_COMPLEX")

        self.assertEqual(hla["evidence_class"], "measured-protein")
        self.assertEqual(hla["measured_dataset_count"], 2)
        self.assertTrue(hla["protein_holdout_eligible"])
        self.assertTrue(hla["donor_holdout_eligible"])
        self.assertFalse(hla["hgsoc_inference_authorized"])

    def test_unmeasured_cd74_is_marked_unvalidated_extrapolation(self):
        rows = self.mod.build_target_coverage_rows(
            measured_targets={"reference": {"HLA_DR_COMPLEX"}},
            protein_entries={"P04233"},
            gene_names={"CD74"},
        )
        cd74 = next(row for row in rows if row["target_id"] == "CD74")

        self.assertEqual(cd74["evidence_class"], "unvalidated-extrapolation")
        self.assertFalse(cd74["protein_holdout_eligible"])
        self.assertFalse(cd74["hgsoc_inference_authorized"])

    def test_measured_target_without_donor_labels_is_not_donor_holdout_eligible(self):
        rows = self.mod.build_target_coverage_rows(
            measured_targets={"GSE200417_CITE": {"CD74"}},
            protein_entries={"P04233"},
            gene_names={"CD74"},
            dataset_donor_counts={"GSE200417_CITE": 0},
        )
        cd74 = next(row for row in rows if row["target_id"] == "CD74")

        self.assertTrue(cd74["protein_holdout_eligible"])
        self.assertFalse(cd74["donor_holdout_eligible"])

    def test_task5_gate_requires_measured_hla_ii_and_complete_embeddings(self):
        self.assertTrue(
            self.mod.task5_gate_passes(
                measured_targets={"ref": {"HLA_DR_COMPLEX"}},
                protein_entries={"P01903", "P01911"},
                gene_names={"HLA-DRA", "HLA-DRB1"},
            )
        )
        self.assertFalse(
            self.mod.task5_gate_passes(
                measured_targets={"ref": {"HLA_DR_COMPLEX"}},
                protein_entries={"P01903"},
                gene_names={"HLA-DRA", "HLA-DRB1"},
            )
        )


if __name__ == "__main__":
    unittest.main()
