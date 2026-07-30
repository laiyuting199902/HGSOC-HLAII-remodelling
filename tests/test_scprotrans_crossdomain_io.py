#!/usr/bin/env python3
"""Tests for cross-domain scProTrans training and inference contracts."""

import importlib.util
import gzip
import hashlib
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import torch


HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT_PATH = os.path.normpath(
    os.path.join(HERE, "..", "scripts", "39_scprotrans_crossdomain_training.py")
)


def load_module():
    spec = importlib.util.spec_from_file_location("scprotrans_crossdomain", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestDataSplits(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_module()

    def test_donor_holdout_has_no_donor_leakage(self):
        donors = np.array(["P1", "P1", "P2", "P2", "P3", "P3"])
        train_idx, test_idx = self.mod.make_donor_holdout(donors, {"P2"})
        self.assertEqual(set(donors[train_idx]), {"P1", "P3"})
        self.assertEqual(set(donors[test_idx]), {"P2"})
        self.assertTrue(set(train_idx).isdisjoint(set(test_idx)))

    def test_empty_or_total_donor_holdout_is_rejected(self):
        donors = np.array(["P1", "P1", "P2", "P2"])
        with self.assertRaises(ValueError):
            self.mod.make_donor_holdout(donors, set())
        with self.assertRaises(ValueError):
            self.mod.make_donor_holdout(donors, {"P1", "P2"})

    def test_protein_holdout_removes_target_from_training(self):
        proteins = ["HLA_DR_COMPLEX", "PD_L1", "CD47"]
        train_idx, test_idx = self.mod.make_protein_holdout(
            proteins, {"HLA_DR_COMPLEX"}
        )
        self.assertEqual([proteins[i] for i in train_idx], ["PD_L1", "CD47"])
        self.assertEqual([proteins[i] for i in test_idx], ["HLA_DR_COMPLEX"])

    def test_group_folds_keep_each_donor_in_exactly_one_test_fold(self):
        donors = np.array(["P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4"])
        folds = self.mod.make_group_folds(donors, n_splits=2, seed=7)
        observed_test = []
        for train_idx, test_idx in folds:
            self.assertTrue(set(donors[train_idx]).isdisjoint(set(donors[test_idx])))
            observed_test.extend(donors[test_idx])
        self.assertCountEqual(observed_test, donors.tolist())


class TestEmbeddingAndModelContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_module()

    def test_complex_embedding_is_mean_of_all_chains(self):
        embeddings = {
            "P01903": np.array([1.0, 3.0, 5.0]),
            "P01911": np.array([3.0, 5.0, 7.0]),
        }
        observed = self.mod.mean_chain_embedding(
            embeddings, ["P01903", "P01911"]
        )
        np.testing.assert_allclose(observed, [2.0, 4.0, 6.0])

    def test_feature_alignment_uses_frozen_order(self):
        matrix = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        observed = self.mod.align_feature_matrix(
            matrix,
            source_names=["B", "A", "C"],
            target_order=["A", "B"],
        )
        np.testing.assert_array_equal(observed, [[2.0, 1.0], [5.0, 4.0]])
        with self.assertRaises(KeyError):
            self.mod.align_feature_matrix(
                matrix,
                source_names=["B", "A", "C"],
                target_order=["A", "MISSING"],
            )

    def test_configurable_protrans_forward_shapes(self):
        model = self.mod.ProTrans(
            protein_dim=3,
            gene_dim=4,
            cell_dim=2,
            hidden_dim=4,
            heads=2,
            head_dim=2,
        )
        protein_embedding = torch.randn(2, 3, 3)
        gene_embedding = torch.randn(2, 5, 4)
        cell_embedding = torch.randn(2, 2)
        expression = torch.rand(2, 5)
        prediction, attention = model(
            protein_embedding,
            gene_embedding,
            cell_embedding,
            expression,
        )
        self.assertEqual(tuple(prediction.shape), (2, 3))
        self.assertEqual(tuple(attention.shape), (2, 3, 5))

    def test_forward_accepts_shared_two_dimensional_sequence_embeddings(self):
        model = self.mod.ProTrans(
            protein_dim=3,
            gene_dim=4,
            cell_dim=2,
            hidden_dim=4,
            heads=2,
            head_dim=2,
        )
        prediction, attention = model(
            torch.randn(3, 3),
            torch.randn(5, 4),
            torch.randn(2, 2),
            torch.rand(2, 5),
        )
        self.assertEqual(tuple(prediction.shape), (2, 3))
        self.assertEqual(tuple(attention.shape), (2, 3, 5))

    def test_factorized_sequence_cell_projection_matches_explicit_linear(self):
        torch.manual_seed(4)
        layer = torch.nn.Linear(7, 5)
        sequence = torch.randn(3, 4)
        cells = torch.randn(2, 3)
        explicit = layer(
            torch.cat(
                [
                    sequence.unsqueeze(0).expand(2, -1, -1),
                    cells[:, None, :].expand(-1, 3, -1),
                ],
                dim=2,
            )
        )
        observed = self.mod.factorized_sequence_cell_projection(
            layer,
            sequence,
            cells,
        )
        torch.testing.assert_close(observed, explicit)

    def test_model_bundle_roundtrip_freezes_orders_and_seed(self):
        model = self.mod.ProTrans(
            protein_dim=3,
            gene_dim=4,
            cell_dim=2,
            hidden_dim=4,
            heads=2,
            head_dim=2,
        )
        metadata = {
            "gene_order": ["A", "B"],
            "protein_order": ["HLA_DR_COMPLEX"],
            "normalization": "library_size_1e4_log1p",
            "seed": 17,
            "model_config": model.config,
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "bundle.pt")
            self.mod.save_model_bundle(path, model, metadata)
            restored_model, restored_metadata = self.mod.load_model_bundle(path)

        self.assertEqual(restored_metadata, metadata)
        for left, right in zip(model.parameters(), restored_model.parameters()):
            self.assertTrue(torch.equal(left, right))


class TestReferenceAndTargetPreprocessing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_module()

    def test_binary_csc_loader_reconstructs_feature_by_cell_matrix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "counts"
            np.array([0, 2, 1], dtype="<u4").tofile(str(prefix) + "_i.bin")
            np.array([1, 2, 3], dtype="<u4").tofile(str(prefix) + "_x.bin")
            (Path(str(prefix) + "_p.tsv")).write_text(
                "column_pointer\n0\n2\n3\n", encoding="ascii"
            )
            (Path(str(prefix) + "_manifest.tsv")).write_text(
                "key\tvalue\nfeatures\t3\ncells\t2\nnonzero\t3\n"
                "endian\tlittle\nindex_base\tzero\nvalue_type\tuint32\n",
                encoding="ascii",
            )
            observed = self.mod.read_binary_csc(prefix)

        np.testing.assert_array_equal(
            observed.toarray(), np.array([[1, 0], [0, 3], [2, 0]])
        )

    def test_library_size_log1p_normalization_preserves_zeroes(self):
        counts = np.array([[1.0, 1.0], [0.0, 4.0]])
        observed = self.mod.normalize_total_log1p(counts, target_sum=10.0)
        expected = np.log1p(np.array([[5.0, 5.0], [0.0, 10.0]]))
        np.testing.assert_allclose(observed, expected, rtol=1e-6)
        self.assertEqual(observed[1, 0], 0.0)

    def test_targetwise_adt_transform_does_not_depend_on_panel_total(self):
        targets = np.array([[0.0, 3.0], [8.0, 15.0]], dtype=float)
        observed = self.mod.transform_target_adt_counts(targets)
        np.testing.assert_allclose(observed, np.log1p(targets))
        with self.assertRaises(ValueError):
            self.mod.transform_target_adt_counts(np.array([[1.0, -1.0]]))

    def test_duplicate_adt_labels_are_averaged_within_canonical_target(self):
        adt = np.array([[2.0, 4.0, 7.0], [6.0, 10.0, 9.0]])
        labels = ["HLA-DR-1", "HLA-DR-2", "CD47"]
        mapper = lambda label: "HLA_DR_COMPLEX" if label.startswith("HLA-DR") else label
        observed, targets = self.mod.aggregate_target_adt(adt, labels, mapper)
        self.assertEqual(targets, ["HLA_DR_COMPLEX", "CD47"])
        np.testing.assert_allclose(observed, [[3.0, 7.0], [8.0, 9.0]])

    def test_gene_selection_keeps_forced_targets_before_variance_rank(self):
        observed = self.mod.select_model_genes(
            reference_genes=["A", "B", "HLA-DRA", "CD74"],
            target_genes=["A", "B", "HLA-DRA", "CD74"],
            embedding_genes={"A", "B", "HLA-DRA", "CD74"},
            variance_scores={"A": 9.0, "B": 8.0, "HLA-DRA": 0.1, "CD74": 0.2},
            forced_genes=["HLA-DRA", "CD74"],
            max_genes=3,
        )
        self.assertEqual(observed, ["HLA-DRA", "CD74", "A"])

    def test_streamed_feature_by_cell_tsv_returns_cells_by_requested_features(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "matrix.tsv.gz"
            with gzip.open(path, "wt", encoding="ascii") as handle:
                handle.write("cell1\tcell2\n")
                handle.write("A\t1\t2\n")
                handle.write("B\t3\t4\n")
                handle.write("C\t5\t6\n")
            cells, observed = self.mod.read_feature_by_cell_tsv(
                path, requested_features=["C", "A"]
            )
        self.assertEqual(cells, ["cell1", "cell2"])
        np.testing.assert_array_equal(observed, [[5, 1], [6, 2]])

    def test_streamed_tsv_skips_feature_barcode_annotation_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "matrix.tsv.gz"
            with gzip.open(path, "wt", encoding="ascii") as handle:
                handle.write("\tbarcode\tcell1\tcell2\n")
                handle.write("HLA-DR\tAATAGC\t6\t18\n")
                handle.write("CD274_PD-L1\tGTTGTC\t4\t0\n")
            cells, observed = self.mod.read_feature_by_cell_tsv(
                path,
                requested_features=["HLA-DR", "CD274_PD-L1"],
            )
        self.assertEqual(cells, ["cell1", "cell2"])
        np.testing.assert_array_equal(observed, [[6, 4], [18, 0]])

    def test_official_gene_embedding_rows_follow_requested_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "gene.npz"
            np.savez(
                path,
                gene=np.array(["A", "B", "C"]),
                embedding=np.array([[1, 2], [3, 4], [5, 6]], dtype=np.float32),
            )
            observed = self.mod.load_gene_embeddings(path, ["C", "A"])
        np.testing.assert_array_equal(observed, [[5, 6], [1, 2]])

    def test_complex_protein_embedding_averages_requested_h5_entries(self):
        import h5py

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "protein.h5"
            with h5py.File(path, "w") as handle:
                handle.create_dataset("P1", data=np.array([1, 3, 5], dtype=np.float16))
                handle.create_dataset("P2", data=np.array([3, 5, 7], dtype=np.float16))
            observed = self.mod.load_protein_embeddings(
                path,
                {"HLA_DR_COMPLEX": ["P1", "P2"]},
            )
        np.testing.assert_allclose(observed, [[2, 4, 6]])


class TestProductionPipelineContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_module()

    def test_h5ad_csr_subset_reads_only_requested_rows_and_columns(self):
        import h5py
        from scipy import sparse

        matrix = sparse.csr_matrix(
            np.array(
                [
                    [1, 0, 2, 0],
                    [0, 3, 0, 4],
                    [5, 0, 6, 0],
                ],
                dtype=np.float32,
            )
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "mini.h5ad"
            with h5py.File(path, "w") as handle:
                layers = handle.create_group("layers")
                counts = layers.create_group("counts")
                counts.create_dataset("data", data=matrix.data)
                counts.create_dataset("indices", data=matrix.indices)
                counts.create_dataset("indptr", data=matrix.indptr)
                counts.attrs["shape"] = matrix.shape
            observed = self.mod.read_h5ad_csr_subset(
                path,
                row_indices=[2, 0],
                column_indices=[2, 0],
                shape=matrix.shape,
            )

        np.testing.assert_array_equal(observed.toarray(), [[6, 5], [2, 1]])

    def test_group_balanced_sampling_is_deterministic_and_capped(self):
        groups = np.array(["A"] * 6 + ["B"] * 3 + ["C"] * 5)
        first = self.mod.sample_indices_by_group(groups, max_per_group=4, seed=17)
        second = self.mod.sample_indices_by_group(groups, max_per_group=4, seed=17)
        np.testing.assert_array_equal(first, second)
        sampled = groups[first]
        self.assertEqual(np.sum(sampled == "A"), 4)
        self.assertEqual(np.sum(sampled == "B"), 3)
        self.assertEqual(np.sum(sampled == "C"), 4)

    def test_sparse_feature_variance_matches_dense_population_variance(self):
        from scipy import sparse

        matrix = np.array([[0, 1, 2], [2, 1, 4], [4, 1, 6]], dtype=np.float32)
        observed = self.mod.sparse_feature_variance(sparse.csr_matrix(matrix))
        np.testing.assert_allclose(observed, np.var(matrix, axis=0), atol=1e-7)

    def test_prepared_block_roundtrip_preserves_sparse_counts_and_metadata(self):
        from scipy import sparse

        block = {
            "counts": sparse.csr_matrix([[1, 0, 2], [0, 3, 0]]),
            "expression": sparse.csr_matrix(
                np.array([[1.0, 0.0], [0.0, 2.0]], dtype=np.float32)
            ),
            "cells": np.array(["c1", "c2"], dtype=object),
            "donors": np.array(["D1", "D2"]),
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest = self.mod.save_prepared_block(Path(tmpdir), "reference", block)
            observed = self.mod.load_prepared_block(Path(tmpdir), "reference")
            for filename, digest in manifest["sha256"].items():
                with open(Path(tmpdir) / filename, "rb") as handle:
                    self.assertEqual(hashlib.sha256(handle.read()).hexdigest(), digest)
        np.testing.assert_array_equal(observed["counts"].toarray(), block["counts"].toarray())
        np.testing.assert_array_equal(observed["expression"], block["expression"].toarray())
        np.testing.assert_array_equal(observed["cells"], block["cells"])
        np.testing.assert_array_equal(observed["donors"], block["donors"])

    def test_patient_paired_effect_uses_patient_not_cell_as_unit(self):
        prediction = np.array([1, 1, 3, 3, 2, 2, 5, 5], dtype=float)
        patients = np.array(["P1"] * 4 + ["P2"] * 4)
        stages = np.array(["naive", "naive", "IDS", "IDS"] * 2)
        observed = self.mod.compute_patient_paired_effect(
            prediction,
            patients,
            stages,
            pre_label="naive",
            post_label="IDS",
        )
        self.assertEqual(observed["patient_count"], 2)
        np.testing.assert_allclose(observed["patient_effects"], [2.0, 3.0])
        self.assertEqual(observed["mean_effect"], 2.5)

    def test_eligible_patient_reader_uses_predefined_minimum_eoc_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "patients.tsv"
            path.write_text(
                "patient_id\teoc_pair_min10\nP1\tTRUE\nP2\tFALSE\nP3\tTRUE\n",
                encoding="ascii",
            )
            observed = self.mod.read_eligible_patients(path)
        self.assertEqual(observed, {"P1", "P3"})

    def test_noninferiority_margin_is_explicit(self):
        self.assertTrue(
            self.mod.not_worse_than_baseline(
                model_metric=0.31, baseline_metric=0.35, margin=0.05
            )
        )
        self.assertFalse(
            self.mod.not_worse_than_baseline(
                model_metric=0.20, baseline_metric=0.35, margin=0.05
            )
        )

    def test_cli_exposes_pilot_and_full_runs(self):
        args = self.mod.parse_args(
            [
                "--run",
                "pilot",
                "--max-genes",
                "256",
                "--seeds",
                "11",
                "--model-hidden-dim",
                "500",
                "--model-heads",
                "10",
            ]
        )
        self.assertEqual(args.run, "pilot")
        self.assertEqual(args.max_genes, 256)
        self.assertEqual(args.seeds, [11])
        self.assertEqual(args.model_hidden_dim, 500)
        self.assertEqual(args.model_heads, 10)

    def test_training_target_mapper_preserves_complex_semantics(self):
        expected = {
            "HLA-DR": "HLA_DR_COMPLEX",
            "HLA.DR": "HLA_DR_COMPLEX",
            "HLA-A-B-C": "HLA_I_PAN",
            "CD274_PD-L1": "PD_L1",
            "CD47": "CD47",
            "Mouse_IgG1": None,
        }
        self.assertEqual(
            {label: self.mod.canonical_training_target(label) for label in expected},
            expected,
        )

    def test_representation_manifest_requires_reference_only_inductive_fit(self):
        manifest = {
            "mode": "reference_only_inductive",
            "fit_donors": ["D1", "D2"],
            "mapped_donors": ["D3"],
            "query_adaptation_epochs": 0,
            "gene_order_sha256": "abc",
        }
        self.mod.validate_representation_manifest(manifest)
        with self.assertRaises(ValueError):
            self.mod.validate_representation_manifest(
                {**manifest, "mode": "joint_reference_target"}
            )
        with self.assertRaises(ValueError):
            self.mod.validate_representation_manifest(
                {**manifest, "fit_donors": ["D1", "D3"]}
            )

    def test_reference_scvi_maps_query_without_query_training(self):
        try:
            import scvi  # noqa: F401
        except ImportError:
            self.skipTest("scvi-tools is not installed")
        rng = np.random.default_rng(21)
        reference = rng.poisson(2.0, size=(24, 6)).astype(np.float32)
        query = rng.poisson(2.0, size=(8, 6)).astype(np.float32)
        with tempfile.TemporaryDirectory() as tmpdir:
            reference_dir = Path(tmpdir) / "reference"
            query_dir = Path(tmpdir) / "query"
            reference_latent = self.mod.fit_reference_scvi_latent(
                reference,
                gene_order=[f"G{i}" for i in range(6)],
                model_dir=reference_dir,
                n_latent=3,
                max_epochs=1,
                seed=7,
                progress=False,
            )
            query_latent = self.mod.map_scvi_query_latent(
                query,
                batches=["query"] * len(query),
                gene_order=[f"G{i}" for i in range(6)],
                reference_model_dir=reference_dir,
                output_model_dir=query_dir,
                max_epochs=0,
                seed=7,
                progress=False,
            )
        self.assertEqual(reference_latent.shape, (24, 3))
        self.assertEqual(query_latent.shape, (8, 3))

    def test_double_gzip_rds_selected_export_preserves_frozen_gene_order(self):
        if shutil.which("Rscript") is None:
            self.skipTest("Rscript is required for the RDS export contract")
        with tempfile.TemporaryDirectory() as tmpdir:
            inner = Path(tmpdir) / "inner.rds.gz"
            outer = Path(tmpdir) / "matrix.rds.gz"
            output = Path(tmpdir) / "selected.tsv.gz"
            r_code = (
                "x <- matrix(1:9, nrow=3, byrow=TRUE, "
                "dimnames=list(c('A','B','C'), c('c1','c2','c3'))); "
                "con <- gzfile(commandArgs(TRUE)[1], 'wb'); saveRDS(x, con); close(con)"
            )
            subprocess.run(
                ["Rscript", "-e", r_code, str(inner)],
                check=True,
                capture_output=True,
                text=True,
            )
            with open(outer, "wb") as handle:
                subprocess.run(["gzip", "-c", inner], check=True, stdout=handle)

            self.mod.export_double_gzip_rds_selected(
                outer,
                gene_order=["C", "A"],
                output_path=output,
            )
            cells, matrix = self.mod.read_feature_by_cell_tsv(
                output,
                requested_features=["C", "A"],
            )

        self.assertEqual(cells, ["c1", "c2", "c3"])
        np.testing.assert_array_equal(matrix, [[7, 1], [8, 2], [9, 3]])


class TestTrainingAndBaselines(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_module()

    def test_training_and_prediction_are_reproducible_for_fixed_seed(self):
        rng = np.random.default_rng(5)
        expression = rng.random((12, 4), dtype=np.float32)
        protein = (expression[:, :2].sum(axis=1, keepdims=True) + 0.2).astype(np.float32)
        latent = rng.normal(size=(12, 2)).astype(np.float32)
        gene_embedding = rng.normal(size=(4, 3)).astype(np.float32)
        protein_embedding = rng.normal(size=(1, 2)).astype(np.float32)
        config = dict(
            protein_dim=2,
            gene_dim=3,
            cell_dim=2,
            hidden_dim=4,
            heads=2,
            head_dim=2,
        )
        predictions = []
        for _ in range(2):
            model, history = self.mod.train_protrans(
                expression,
                protein,
                latent,
                gene_embedding,
                protein_embedding,
                model_config=config,
                seed=13,
                epochs=3,
                batch_size=4,
                learning_rate=0.01,
                patience=3,
                device="cpu",
            )
            self.assertGreaterEqual(len(history), 1)
            predictions.append(
                self.mod.predict_protrans(
                    model,
                    expression,
                    latent,
                    gene_embedding,
                    protein_embedding,
                    batch_size=5,
                    device="cpu",
                )
            )
        np.testing.assert_allclose(predictions[0], predictions[1], atol=1e-7)

    def test_output_bias_initialization_prevents_dead_final_relu(self):
        model = self.mod.ProTrans(
            protein_dim=2,
            gene_dim=3,
            cell_dim=2,
            hidden_dim=4,
            heads=2,
            head_dim=2,
        )
        protein_truth = np.array([[0.0], [2.0], [4.0]], dtype=np.float32)
        self.mod.initialize_output_bias(model, protein_truth)
        self.assertGreater(float(model.linear_2.bias.detach()), 0.0)
        self.assertAlmostEqual(float(model.linear_2.bias.detach()), 2.0)

    def test_ridge_baseline_recovers_a_linear_signal(self):
        x = np.arange(30, dtype=float).reshape(10, 3)
        y = (2 * x[:, [0]] - x[:, [1]])
        model = self.mod.fit_ridge_baseline(x, y, alpha=1e-8)
        prediction = self.mod.predict_ridge_baseline(model, x)
        self.assertGreater(np.corrcoef(y[:, 0], prediction[:, 0])[0, 1], 0.999)

    def test_cognate_rna_baseline_averages_frozen_target_genes(self):
        expression = np.array([[1.0, 3.0, 9.0], [2.0, 6.0, 8.0]])
        observed = self.mod.cognate_rna_score(
            expression,
            gene_order=["HLA-DRA", "HLA-DRB1", "OTHER"],
            target="HLA_DR_COMPLEX",
        )
        np.testing.assert_allclose(observed, [2.0, 4.0])

    def test_group_metric_records_keep_each_donor_separate(self):
        truth = np.arange(12, dtype=float)
        prediction = truth.copy()
        groups = np.array(["D1"] * 4 + ["D2"] * 4 + ["D3"] * 4)
        records = self.mod.compute_group_metric_records(
            truth,
            prediction,
            groups,
            metric="spearman",
        )
        self.assertEqual([record["group"] for record in records], ["D1", "D2", "D3"])
        self.assertTrue(all(record["value"] == 1.0 for record in records))


class TestMetricsOODAndGate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_module()

    def test_prediction_metrics_include_required_statistics(self):
        truth = np.array([0.0, 1.0, 2.0, 3.0])
        prediction = np.array([0.1, 1.1, 1.9, 3.1])
        metrics = self.mod.compute_prediction_metrics(truth, prediction)
        self.assertEqual(
            set(metrics),
            {
                "spearman",
                "pearson",
                "mae",
                "cosine",
                "calibration_intercept",
                "calibration_slope",
            },
        )
        self.assertGreater(metrics["spearman"], 0.9)
        self.assertLess(metrics["mae"], 0.2)

    def test_calibration_reports_intercept_instead_of_forcing_origin(self):
        prediction = np.arange(6, dtype=float)
        truth = 3.0 + 2.0 * prediction
        metrics = self.mod.compute_prediction_metrics(truth, prediction)
        self.assertAlmostEqual(metrics["calibration_intercept"], 3.0)
        self.assertAlmostEqual(metrics["calibration_slope"], 2.0)

    def test_ood_model_marks_far_target_as_out_of_domain(self):
        reference = np.array(
            [[0.0, 0.0], [0.1, 0.0], [0.0, 0.1], [0.1, 0.1], [0.05, 0.05]]
        )
        model = self.mod.fit_knn_ood(reference, k=2, quantile=0.95)
        scores, in_domain = self.mod.score_knn_ood(
            model,
            np.array([[0.05, 0.05], [10.0, 10.0]]),
        )
        self.assertTrue(in_domain[0])
        self.assertFalse(in_domain[1])
        self.assertLess(scores[0], scores[1])

    def test_ood_threshold_uses_leave_one_donor_out_neighbors(self):
        reference = np.array(
            [
                [0.00, 0.00],
                [0.01, 0.00],
                [1.00, 0.00],
                [1.01, 0.00],
                [2.00, 0.00],
                [2.01, 0.00],
            ]
        )
        donors = np.array(["D1", "D1", "D2", "D2", "D3", "D3"])
        model = self.mod.fit_knn_ood(
            reference,
            reference_groups=donors,
            k=1,
            quantile=0.95,
        )
        self.assertGreater(model["threshold"], 0.9)
        self.assertEqual(model["calibration"], "leave_one_group_out")

    def test_patient_equal_ood_fraction_is_not_cell_weighted(self):
        in_domain = np.array([True] * 90 + [False] * 10 + [False] * 10)
        patients = np.array(["P1"] * 100 + ["P2"] * 10)
        observed = self.mod.patient_equal_in_domain_fraction(in_domain, patients)
        self.assertAlmostEqual(observed["cell_fraction"], 90 / 110)
        self.assertAlmostEqual(observed["patient_equal_fraction"], 0.45)

    def test_gate_requires_all_predefined_conditions(self):
        passing = self.mod.evaluate_gate3(
            donor_summary={
                "estimate": 0.45,
                "ci_low": 0.12,
                "donor_count": 9,
                "finite_donor_count": 9,
                "finite_draws": 980,
                "iterations": 1000,
            },
            independent_model_metrics=[0.40, 0.42, 0.38, 0.41],
            independent_baseline_metrics=[0.39, 0.40, 0.37, 0.39],
            hgsoc_in_domain=[True] * 8 + [False] * 2,
            hgsoc_patients=["P1"] * 5 + ["P2"] * 5,
            seed_patient_effects=np.array(
                [
                    [0.2, 0.4, 0.6],
                    [0.3, 0.5, 0.7],
                    [0.1, 0.4, 0.8],
                    [0.2, 0.5, 0.7],
                    [0.3, 0.4, 0.6],
                ]
            ),
            bundle_audit={
                "roundtrip_match": True,
                "artifact_hashes_match": True,
                "encoder_present": True,
                "split_manifest_present": True,
            },
            bootstrap_iterations=200,
            seed=17,
        )
        self.assertTrue(passing["gate3_passed"])
        self.assertEqual(passing["failed_conditions"], [])

        failing = self.mod.evaluate_gate3(
            donor_summary={
                "estimate": 0.45,
                "ci_low": 0.12,
                "donor_count": 9,
                "finite_donor_count": 9,
                "finite_draws": 980,
                "iterations": 1000,
            },
            independent_model_metrics=[0.40, 0.42, 0.38, 0.41],
            independent_baseline_metrics=[0.39, 0.40, 0.37, 0.39],
            hgsoc_in_domain=[True] * 9 + [False] * 11,
            hgsoc_patients=["P1"] * 10 + ["P2"] * 10,
            seed_patient_effects=np.array(
                [
                    [0.2, 0.4, 0.6],
                    [0.3, 0.5, 0.7],
                    [0.1, 0.4, 0.8],
                    [0.2, 0.5, 0.7],
                    [0.3, 0.4, 0.6],
                ]
            ),
            bundle_audit={
                "roundtrip_match": True,
                "artifact_hashes_match": True,
                "encoder_present": True,
                "split_manifest_present": True,
            },
            bootstrap_iterations=200,
            seed=17,
        )
        self.assertFalse(failing["gate3_passed"])
        self.assertIn("hgsoc_ood_coverage", failing["failed_conditions"])

    def test_gate_rejects_invalid_summary_ranges_and_unverified_bundle(self):
        with self.assertRaises(ValueError):
            self.mod.evaluate_gate3(
                donor_summary={
                    "estimate": 0.45,
                    "ci_low": 0.12,
                    "donor_count": 9,
                    "finite_donor_count": 9,
                    "finite_draws": 1001,
                    "iterations": 1000,
                },
                independent_model_metrics=[0.4, 0.4, 0.4],
                independent_baseline_metrics=[0.4, 0.4, 0.4],
                hgsoc_in_domain=[True, True],
                hgsoc_patients=["P1", "P2"],
                seed_patient_effects=np.ones((5, 2)),
                bundle_audit={},
            )

    def test_gate_rejects_excess_undefined_donor_metrics(self):
        result = self.mod.evaluate_gate3(
            donor_summary={
                "estimate": 0.62,
                "ci_low": 0.53,
                "donor_count": 9,
                "finite_donor_count": 5,
                "finite_draws": 1000,
                "iterations": 1000,
            },
            independent_model_metrics=[0.4, 0.4, 0.4],
            independent_baseline_metrics=[0.4, 0.4, 0.4],
            hgsoc_in_domain=[True, True],
            hgsoc_patients=["P1", "P2"],
            seed_patient_effects=np.array(
                [[0.1, 0.3], [0.1, 0.3], [0.1, 0.3], [0.1, 0.3], [0.1, 0.3]]
            ),
            bundle_audit={
                "roundtrip_match": True,
                "artifact_hashes_match": True,
                "encoder_present": True,
                "split_manifest_present": True,
            },
            bootstrap_iterations=100,
        )
        self.assertFalse(result["donor_holdout_performance"])

    def test_paired_noninferiority_uses_ci_not_only_point_estimate(self):
        observed = self.mod.paired_noninferiority_summary(
            model_metrics=[0.40, 0.40, 0.40, 0.10],
            baseline_metrics=[0.40, 0.40, 0.40, 0.20],
            margin=0.05,
            iterations=1000,
            seed=9,
        )
        self.assertGreater(observed["mean_difference"], -0.05)
        self.assertLessEqual(observed["ci_low"], -0.05)
        self.assertFalse(observed["passed"])

    def test_seed_stability_is_computed_from_patient_by_seed_matrix(self):
        observed = self.mod.summarize_seed_patient_effects(
            np.array(
                [
                    [0.1, 0.5, 0.9],
                    [0.2, 0.5, 0.8],
                    [0.1, 0.4, 0.9],
                    [0.2, 0.4, 0.8],
                    [0.1, 0.5, 0.8],
                ]
            )
        )
        self.assertEqual(observed["seed_count"], 5)
        self.assertTrue(observed["all_seed_means_positive"])
        self.assertLess(observed["seed_effect_sd"], observed["patient_effect_sd"])

    def test_critical_gate_failure_stops_before_five_seed_run(self):
        result = self.mod.evaluate_gate3_critical_stage(
            donor_summary={
                "estimate": 0.63,
                "ci_low": 0.57,
                "donor_count": 9,
                "finite_donor_count": 9,
                "finite_draws": 1000,
                "iterations": 1000,
            },
            independent_model_metrics=[0.47, 0.31, -0.03],
            independent_baseline_metrics=[0.76, 0.36, 0.24],
            hgsoc_patient_equal_in_domain_fraction=0.148,
            bundle_audit={
                "roundtrip_match": True,
                "artifact_hashes_match": True,
                "encoder_present": True,
                "split_manifest_present": True,
            },
            bootstrap_iterations=200,
            seed=5,
        )
        self.assertFalse(result["proceed_to_five_seed"])
        self.assertIn("independent_transfer", result["failed_conditions"])
        self.assertIn("hgsoc_ood_coverage", result["failed_conditions"])

    def test_cluster_bootstrap_resamples_donors_not_individual_cells(self):
        truth = np.arange(8, dtype=float)
        prediction = truth.copy()
        donors = np.array(["P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4"])
        observed = self.mod.cluster_bootstrap_metric(
            truth,
            prediction,
            donors,
            metric="spearman",
            iterations=100,
            seed=11,
        )
        self.assertAlmostEqual(observed["estimate"], 1.0)
        self.assertAlmostEqual(observed["ci_low"], 1.0)
        self.assertAlmostEqual(observed["ci_high"], 1.0)

    def test_cluster_bootstrap_estimate_weights_donors_equally(self):
        truth = np.concatenate(
            [np.arange(100, dtype=float), np.arange(5), np.arange(5), np.arange(5)]
        )
        prediction = np.concatenate(
            [np.arange(100, dtype=float)[::-1], np.arange(5), np.arange(5), np.arange(5)]
        )
        donors = np.array(["D1"] * 100 + ["D2"] * 5 + ["D3"] * 5 + ["D4"] * 5)
        observed = self.mod.cluster_bootstrap_metric(
            truth,
            prediction,
            donors,
            metric="spearman",
            iterations=200,
            seed=3,
        )
        self.assertAlmostEqual(observed["estimate"], 0.5)
        self.assertEqual(observed["donor_count"], 4)


if __name__ == "__main__":
    unittest.main()
