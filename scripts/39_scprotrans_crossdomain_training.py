#!/usr/bin/env python3
"""Cross-domain scProTrans contracts, model, metrics, and Gate 3 helpers."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import random
import struct
import subprocess
import zipfile
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy import sparse
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.neighbors import NearestNeighbors


TARGET_ORDER = ("HLA_DR_COMPLEX", "HLA_I_PAN", "PD_L1", "CD47")
TARGET_PROTEIN_ENTRIES = {
    "HLA_DR_COMPLEX": ("P01903", "P01911"),
    "HLA_I_PAN": ("P04439", "P01889", "P10321", "P61769"),
    "PD_L1": ("Q9NZQ7",),
    "CD47": ("Q08722",),
}
TARGET_COGNATE_GENES = {
    "HLA_DR_COMPLEX": ("HLA-DRA", "HLA-DRB1"),
    "HLA_I_PAN": ("HLA-A", "HLA-B", "HLA-C", "B2M"),
    "PD_L1": ("CD274",),
    "CD47": ("CD47",),
}
FORCED_GENES = (
    "HLA-DRA",
    "HLA-DRB1",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "CD74",
    "CIITA",
    "CD274",
    "CD47",
    "EPCAM",
    "PTPRC",
)


def sha256_file(path: os.PathLike[str] | str, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def save_prepared_block(
    directory: os.PathLike[str] | str,
    name: str,
    block: Mapping[str, object],
) -> dict[str, object]:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    if "counts" not in block:
        raise ValueError("Prepared blocks require a counts matrix.")
    counts_path = directory / f"{name}_counts.npz"
    arrays_path = directory / f"{name}_arrays.npz"
    sparse.save_npz(counts_path, sparse.csr_matrix(block["counts"]), compressed=True)
    arrays = {}
    for key, value in block.items():
        if key == "counts":
            continue
        array = value.toarray() if sparse.issparse(value) else np.asarray(value)
        if array.dtype == object:
            array = array.astype(str)
        arrays[key] = array
    np.savez_compressed(arrays_path, **arrays)
    return {
        "name": name,
        "counts_shape": list(sparse.csr_matrix(block["counts"]).shape),
        "array_keys": sorted(arrays),
        "sha256": {
            counts_path.name: sha256_file(counts_path),
            arrays_path.name: sha256_file(arrays_path),
        },
    }


def load_prepared_block(
    directory: os.PathLike[str] | str,
    name: str,
) -> dict[str, object]:
    directory = Path(directory)
    block: dict[str, object] = {
        "counts": sparse.load_npz(directory / f"{name}_counts.npz").tocsr()
    }
    with np.load(directory / f"{name}_arrays.npz", allow_pickle=False) as archive:
        block.update({key: np.asarray(archive[key]) for key in archive.files})
    return block


def canonical_training_target(label: str) -> str | None:
    compact = "".join(character for character in str(label).upper() if character.isalnum())
    if not compact or "IGG" in compact or "CONTROL" in compact or "CTRL" in compact:
        return None
    if "HLADRDPDQ" in compact:
        return "HLA_II_PAN"
    if "HLAABC" in compact:
        return "HLA_I_PAN"
    if compact.startswith("HLADR") or "MHCII" in compact:
        return "HLA_DR_COMPLEX"
    if "PDL1" in compact or compact.startswith("CD274"):
        return "PD_L1"
    if compact.startswith("CD47"):
        return "CD47"
    return None


def make_donor_holdout(
    donors: Sequence[str],
    held_out_donors: set[str],
) -> tuple[np.ndarray, np.ndarray]:
    donors = np.asarray(donors, dtype=str)
    observed = set(donors)
    if not held_out_donors or not held_out_donors.issubset(observed):
        raise ValueError("Held-out donors must be a non-empty subset of observed donors.")
    if held_out_donors == observed:
        raise ValueError("Donor holdout must retain at least one training donor.")
    test_mask = np.isin(donors, sorted(held_out_donors))
    train_idx = np.flatnonzero(~test_mask)
    test_idx = np.flatnonzero(test_mask)
    if set(donors[train_idx]).intersection(set(donors[test_idx])):
        raise AssertionError("Donor leakage detected.")
    return train_idx, test_idx


def make_protein_holdout(
    proteins: Sequence[str],
    held_out_proteins: set[str],
) -> tuple[np.ndarray, np.ndarray]:
    proteins = list(proteins)
    observed = set(proteins)
    if not held_out_proteins or not held_out_proteins.issubset(observed):
        raise ValueError("Held-out proteins must be a non-empty subset of the panel.")
    if held_out_proteins == observed:
        raise ValueError("Protein holdout must retain at least one training protein.")
    train_idx = np.array(
        [index for index, protein in enumerate(proteins) if protein not in held_out_proteins],
        dtype=int,
    )
    test_idx = np.array(
        [index for index, protein in enumerate(proteins) if protein in held_out_proteins],
        dtype=int,
    )
    return train_idx, test_idx


def make_group_folds(
    groups: Sequence[str],
    n_splits: int,
    seed: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    groups = np.asarray(groups, dtype=str)
    unique_groups = np.unique(groups)
    if n_splits < 2 or n_splits > unique_groups.size:
        raise ValueError("n_splits must be between two and the number of groups.")
    rng = np.random.default_rng(seed)
    shuffled = unique_groups.copy()
    rng.shuffle(shuffled)
    folds = []
    for test_groups in np.array_split(shuffled, n_splits):
        test_mask = np.isin(groups, test_groups)
        train_idx = np.flatnonzero(~test_mask)
        test_idx = np.flatnonzero(test_mask)
        if set(groups[train_idx]).intersection(set(groups[test_idx])):
            raise AssertionError("Group leakage detected.")
        folds.append((train_idx, test_idx))
    return folds


def read_binary_csc(prefix: os.PathLike[str] | str) -> sparse.csc_matrix:
    prefix = os.fspath(prefix)
    manifest_path = prefix + "_manifest.tsv"
    with open(manifest_path, newline="", encoding="ascii") as handle:
        manifest = {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
    required = {
        "features",
        "cells",
        "nonzero",
        "endian",
        "index_base",
        "value_type",
    }
    missing = required.difference(manifest)
    if missing:
        raise ValueError(f"CSC manifest is missing fields: {sorted(missing)}")
    if manifest["endian"] != "little" or manifest["index_base"] != "zero":
        raise ValueError("Only little-endian, zero-based CSC checkpoints are supported.")
    features = int(manifest["features"])
    cells = int(manifest["cells"])
    nonzero = int(manifest["nonzero"])
    indices = np.memmap(prefix + "_i.bin", dtype="<u4", mode="r", shape=(nonzero,))
    values = np.memmap(prefix + "_x.bin", dtype="<u4", mode="r", shape=(nonzero,))
    with open(prefix + "_p.tsv", newline="", encoding="ascii") as handle:
        pointers = np.fromiter(
            (int(row["column_pointer"]) for row in csv.DictReader(handle, delimiter="\t")),
            dtype=np.int64,
            count=cells + 1,
        )
    if pointers.size != cells + 1 or pointers[0] != 0 or pointers[-1] != nonzero:
        raise ValueError("CSC column pointers do not match manifest dimensions.")
    return sparse.csc_matrix(
        (values, indices.astype(np.int32, copy=False), pointers.astype(np.int32, copy=False)),
        shape=(features, cells),
    )


def read_h5ad_csr_subset(
    path: os.PathLike[str] | str,
    *,
    row_indices: Sequence[int],
    column_indices: Sequence[int],
    shape: tuple[int, int],
    layer: str = "counts",
) -> sparse.csr_matrix:
    """Load a sparse H5AD layer once, then retain requested rows and columns."""
    import h5py

    rows = np.asarray(row_indices, dtype=np.int64)
    columns = np.asarray(column_indices, dtype=np.int64)
    if rows.ndim != 1 or columns.ndim != 1 or rows.size == 0 or columns.size == 0:
        raise ValueError("H5AD subset indices must be non-empty one-dimensional arrays.")
    with h5py.File(path, "r") as handle:
        group = handle["X"] if layer == "X" else handle["layers"][layer]
        matrix = sparse.csr_matrix(
            (
                np.asarray(group["data"]),
                np.asarray(group["indices"], dtype=np.int32),
                np.asarray(group["indptr"], dtype=np.int64),
            ),
            shape=shape,
        )
    return matrix[rows][:, columns].tocsr()


def sample_indices_by_group(
    groups: Sequence[str],
    *,
    max_per_group: int,
    seed: int,
) -> np.ndarray:
    groups = np.asarray(groups, dtype=str)
    if groups.ndim != 1 or groups.size == 0 or max_per_group < 1:
        raise ValueError("Groups must be non-empty and max_per_group must be positive.")
    rng = np.random.default_rng(seed)
    sampled = []
    for group in np.unique(groups):
        candidates = np.flatnonzero(groups == group)
        if candidates.size > max_per_group:
            candidates = rng.choice(candidates, size=max_per_group, replace=False)
        sampled.extend(np.asarray(candidates, dtype=int).tolist())
    return np.asarray(sorted(sampled), dtype=int)


def sparse_feature_variance(matrix) -> np.ndarray:
    values = sparse.csr_matrix(matrix, dtype=np.float64)
    mean = np.asarray(values.mean(axis=0)).reshape(-1)
    mean_square = np.asarray(values.multiply(values).mean(axis=0)).reshape(-1)
    return np.maximum(mean_square - mean * mean, 0.0)


def read_feature_by_cell_tsv(
    path: os.PathLike[str] | str,
    *,
    requested_features: Sequence[str],
) -> tuple[list[str], np.ndarray]:
    requested = list(requested_features)
    if not requested or len(requested) != len(set(requested)):
        raise ValueError("Requested features must be non-empty and unique.")
    opener = gzip.open if os.fspath(path).endswith(".gz") else open
    rows: dict[str, np.ndarray] = {}
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if header and header[0] == "":
            header = header[1:]
        has_feature_barcode = bool(header and header[0].strip().lower() == "barcode")
        if has_feature_barcode:
            header = header[1:]
        if not header:
            raise ValueError("Feature-by-cell matrix has an empty header.")
        requested_set = set(requested)
        for line in handle:
            label, separator, values_text = line.rstrip("\r\n").partition("\t")
            if not separator or label not in requested_set:
                continue
            if has_feature_barcode:
                _, separator, values_text = values_text.partition("\t")
                if not separator:
                    raise ValueError(f"Feature barcode annotation is missing values: {label}")
            values = np.fromstring(values_text, sep="\t", dtype=np.float32)
            if values.size != len(header):
                raise ValueError(f"Feature row has the wrong cell count: {label}")
            rows[label] = values
    missing = [feature for feature in requested if feature not in rows]
    if missing:
        raise KeyError(f"Requested matrix features were not found: {missing}")
    return header, np.column_stack([rows[feature] for feature in requested])


def export_double_gzip_rds_selected(
    path: os.PathLike[str] | str,
    *,
    gene_order: Sequence[str],
    output_path: os.PathLike[str] | str,
) -> Path:
    """Export frozen genes from a doubly gzip-wrapped RDS matrix as TSV.gz."""
    genes = list(gene_order)
    if not genes or len(genes) != len(set(genes)):
        raise ValueError("The frozen RDS gene order must be non-empty and unique.")
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    gene_path = output.with_name(output.name + ".genes.txt")
    with open(gene_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(genes) + "\n")
    code = """
args <- commandArgs(TRUE)
con <- gzcon(gzfile(args[[1]], "rb"))
x <- readRDS(con)
close(con)
genes <- readLines(args[[2]], warn=FALSE)
missing <- setdiff(genes, rownames(x))
if (length(missing) > 0) stop(paste("Missing RDS genes:", paste(missing, collapse=",")))
m <- as.matrix(x[genes, , drop=FALSE])
out <- gzfile(args[[3]], "wt")
write.table(m, file=out, sep="\t", quote=FALSE, col.names=NA)
close(out)
"""
    try:
        subprocess.run(
            ["Rscript", "-e", code, os.fspath(path), str(gene_path), str(output)],
            check=True,
            capture_output=True,
            text=True,
        )
    finally:
        gene_path.unlink(missing_ok=True)
    return output


def _open_text(path: os.PathLike[str] | str):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if os.fspath(path).endswith(
        ".gz"
    ) else open(path, "r", encoding="utf-8", newline="")


def read_feature_names(path: os.PathLike[str] | str) -> list[str]:
    names = []
    with _open_text(path) as handle:
        next(handle)
        for line in handle:
            name = line.split("\t", 1)[0].strip()
            if name:
                names.append(name)
    return names


def read_line_names(path: os.PathLike[str] | str) -> list[str]:
    with _open_text(path) as handle:
        return [line.rstrip("\r\n").split("\t", 1)[0] for line in handle if line.strip()]


def read_all_feature_by_cell_tsv(
    path: os.PathLike[str] | str,
) -> tuple[list[str], list[str], np.ndarray]:
    features = read_feature_names(path)
    cells, matrix = read_feature_by_cell_tsv(path, requested_features=features)
    return features, cells, matrix


def read_double_gzip_rds_gene_names(path: os.PathLike[str] | str) -> list[str]:
    code = """
args <- commandArgs(TRUE)
con <- gzcon(gzfile(args[[1]], "rb"))
x <- readRDS(con)
close(con)
writeLines(rownames(x))
"""
    result = subprocess.run(
        ["Rscript", "-e", code, os.fspath(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return [line for line in result.stdout.splitlines() if line]


def read_gene_embedding_names(path: os.PathLike[str] | str) -> set[str]:
    shape, fortran_order, dtype, offset = _stored_npy_layout(path, "gene.npy")
    if fortran_order or len(shape) != 1:
        raise ValueError("Official gene names must be a one-dimensional C-order array.")
    values = np.memmap(path, mode="r", dtype=dtype, offset=offset, shape=shape, order="C")
    return {str(value) for value in values}


def inspect_gse194122(path: os.PathLike[str] | str) -> dict[str, object]:
    import anndata as ad

    adata = ad.read_h5ad(path, backed="r")
    try:
        feature_types = adata.var["feature_types"].astype(str).to_numpy()
        var_names = np.asarray(adata.var_names.astype(str))
        gex_indices = np.flatnonzero(feature_types == "GEX")
        adt_indices = np.flatnonzero(feature_types == "ADT")
        return {
            "shape": (int(adata.n_obs), int(adata.n_vars)),
            "cell_names": np.asarray(adata.obs_names.astype(str)),
            "gene_names": var_names[gex_indices].tolist(),
            "gene_indices": gex_indices,
            "adt_labels": var_names[adt_indices].tolist(),
            "adt_indices": adt_indices,
            "donors": adata.obs["DonorNumber"].astype(str).to_numpy(),
            "batches": adata.obs["batch"].astype(str).to_numpy(),
        }
    finally:
        adata.file.close()


def _reorder_targets(
    values: np.ndarray,
    source_targets: Sequence[str],
    target_order: Sequence[str],
) -> np.ndarray:
    positions = {target: index for index, target in enumerate(source_targets)}
    missing = [target for target in target_order if target not in positions]
    if missing:
        raise KeyError(f"Protein targets are missing from the measured panel: {missing}")
    return values[:, [positions[target] for target in target_order]]


def _align_cells(
    left_cells: Sequence[str],
    left_matrix,
    right_cells: Sequence[str],
    right_matrix,
) -> tuple[list[str], object, object]:
    left_index = {cell: index for index, cell in enumerate(left_cells)}
    right_index = {cell: index for index, cell in enumerate(right_cells)}
    common = [cell for cell in left_cells if cell in right_index]
    if not common:
        raise ValueError("The paired modalities do not share any cell identifiers.")
    left_positions = [left_index[cell] for cell in common]
    right_positions = [right_index[cell] for cell in common]
    return common, left_matrix[left_positions], right_matrix[right_positions]


def normalize_total_log1p(matrix, target_sum: float = 10000.0):
    if target_sum <= 0:
        raise ValueError("target_sum must be positive.")
    if sparse.issparse(matrix):
        normalized = matrix.astype(np.float32, copy=True).tocsr()
        totals = np.asarray(normalized.sum(axis=1)).reshape(-1)
        factors = np.divide(
            target_sum,
            totals,
            out=np.zeros_like(totals, dtype=np.float32),
            where=totals > 0,
        )
        normalized = sparse.diags(factors).dot(normalized).tocsr()
        normalized.data = np.log1p(normalized.data)
        return normalized
    array = np.asarray(matrix, dtype=np.float32)
    totals = array.sum(axis=1, keepdims=True)
    factors = np.divide(
        target_sum,
        totals,
        out=np.zeros_like(totals, dtype=np.float32),
        where=totals > 0,
    )
    return np.log1p(array * factors)


def transform_target_adt_counts(matrix) -> np.ndarray:
    """Apply a target-wise transform without using the rest of an ADT panel."""
    values = np.asarray(matrix, dtype=np.float32)
    if values.ndim != 2 or not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("ADT target counts must be a finite non-negative matrix.")
    return np.log1p(values).astype(np.float32, copy=False)


def aggregate_target_adt(
    adt_matrix,
    labels: Sequence[str],
    target_mapper,
) -> tuple[np.ndarray, list[str]]:
    values = np.asarray(adt_matrix, dtype=np.float32)
    if values.ndim != 2 or values.shape[1] != len(labels):
        raise ValueError("ADT matrix columns must align with labels.")
    target_positions: dict[str, list[int]] = {}
    for index, label in enumerate(labels):
        target = target_mapper(label)
        if target is not None:
            target_positions.setdefault(str(target), []).append(index)
    targets = list(target_positions)
    if not targets:
        raise ValueError("No canonical target ADTs were found.")
    aggregated = np.column_stack(
        [values[:, target_positions[target]].mean(axis=1) for target in targets]
    )
    return aggregated.astype(np.float32, copy=False), targets


def select_model_genes(
    *,
    reference_genes: Sequence[str],
    target_genes: Sequence[str],
    embedding_genes: set[str],
    variance_scores: Mapping[str, float],
    forced_genes: Sequence[str],
    max_genes: int,
) -> list[str]:
    if max_genes < 1:
        raise ValueError("max_genes must be positive.")
    common = set(reference_genes).intersection(target_genes, embedding_genes)
    forced = []
    for gene in forced_genes:
        if gene in common and gene not in forced:
            forced.append(gene)
    ranked = sorted(
        common.difference(forced),
        key=lambda gene: (-float(variance_scores.get(gene, float("-inf"))), gene),
    )
    selected = forced + ranked[: max(0, max_genes - len(forced))]
    if len(selected) > max_genes:
        raise ValueError("The forced gene list exceeds max_genes.")
    return selected


def load_hgsoc_selected_counts(
    hgsoc_root: Path,
    *,
    gene_order: Sequence[str],
    max_cells: int,
    seed: int,
    identity_table: Path,
) -> dict[str, object]:
    feature_path = hgsoc_root / "GSE266577_seurat_features.txt.gz"
    barcode_path = hgsoc_root / "GSE266577_barcodes.txt.gz"
    metadata_path = hgsoc_root / "GSE266577_metadata.txt.gz"
    prefix = hgsoc_root / "derived" / "eoc_csc" / "gse266577_eoc_counts"
    subset_map_path = hgsoc_root / "derived" / "eoc_csc" / "gse266577_eoc_subset_map.tsv"

    features = read_line_names(feature_path)
    feature_index = {gene: index for index, gene in enumerate(features)}
    missing = [gene for gene in gene_order if gene not in feature_index]
    if missing:
        raise KeyError(f"Frozen genes are missing from GSE266577: {missing}")
    full_matrix = read_binary_csc(prefix)
    counts = full_matrix[[feature_index[gene] for gene in gene_order], :].T.tocsr()

    barcodes = read_line_names(barcode_path)
    cell_indices = []
    with open(subset_map_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cell_indices.append(int(row["cell_index"]))
    if len(cell_indices) != counts.shape[0]:
        raise ValueError("GSE266577 EOC subset map does not align with the CSC columns.")
    offset = 1 if min(cell_indices) >= 1 and max(cell_indices) <= len(barcodes) else 0
    cell_names = np.asarray([barcodes[index - offset] for index in cell_indices], dtype=str)

    metadata = {}
    with gzip.open(metadata_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            metadata[row["cell_name"]] = (
                row["publication_patient_code_final"],
                row["publication_sample_code_final"],
                row["treatment_stage"],
            )
    missing_metadata = [cell for cell in cell_names if cell not in metadata]
    if missing_metadata:
        raise ValueError(f"Missing GSE266577 metadata for {len(missing_metadata)} EOC cells.")
    patients = np.asarray([metadata[cell][0] for cell in cell_names], dtype=str)
    samples = np.asarray([metadata[cell][1] for cell in cell_names], dtype=str)
    stages = np.asarray([metadata[cell][2] for cell in cell_names], dtype=str)

    high_confidence_cells: set[str] = set()
    if identity_table.exists():
        with gzip.open(identity_table, "rt", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("strict_eoc") == "TRUE":
                    high_confidence_cells.add(row["cell"])
    high_confidence = np.asarray(
        [cell in high_confidence_cells for cell in cell_names], dtype=bool
    )

    if max_cells > 0 and counts.shape[0] > max_cells:
        per_group = max(1, int(math.ceil(max_cells / np.unique(samples).size)))
        selected = sample_indices_by_group(samples, max_per_group=per_group, seed=seed)
        selected = selected[:max_cells]
        counts = counts[selected]
        cell_names = cell_names[selected]
        patients = patients[selected]
        samples = samples[selected]
        stages = stages[selected]
        high_confidence = high_confidence[selected]
    return {
        "counts": counts,
        "expression": normalize_total_log1p(counts),
        "cells": cell_names,
        "patients": patients,
        "samples": samples,
        "stages": stages,
        "high_confidence": high_confidence,
    }


def load_gse128639_query(
    data_root: Path,
    *,
    gene_order: Sequence[str],
    max_cells: int,
    seed: int,
) -> dict[str, object]:
    root = data_root / "GSE128639"
    rna_path = root / "GSM3681518_MNC_RNA_counts.tsv.gz"
    adt_path = root / "GSM3681519_MNC_ADT_counts.tsv.gz"
    rna_cells, counts = read_feature_by_cell_tsv(rna_path, requested_features=gene_order)
    adt_labels, adt_cells, adt_counts = read_all_feature_by_cell_tsv(adt_path)
    cells, counts, adt_counts = _align_cells(rna_cells, counts, adt_cells, adt_counts)
    if max_cells > 0 and len(cells) > max_cells:
        rng = np.random.default_rng(seed)
        selected = np.sort(rng.choice(len(cells), size=max_cells, replace=False))
        cells = [cells[index] for index in selected]
        counts = counts[selected]
        adt_counts = adt_counts[selected]
    protein, targets = aggregate_target_adt(
        transform_target_adt_counts(adt_counts), adt_labels, canonical_training_target
    )
    target_order = ["HLA_DR_COMPLEX"]
    return {
        "counts": sparse.csr_matrix(counts),
        "expression": normalize_total_log1p(counts),
        "protein": _reorder_targets(protein, targets, target_order),
        "targets": target_order,
        "cells": np.asarray(cells, dtype=str),
        "batches": np.asarray(["GSE128639_BMMC"] * len(cells), dtype=str),
    }


def load_gse254985_query(
    data_root: Path,
    *,
    gene_order: Sequence[str],
    cache_dir: Path,
    max_cells_per_arm: int,
    seed: int,
) -> dict[str, object]:
    root = data_root / "GSE254985" / "extracted"
    configs = [
        (
            "untreated",
            root / "GSM8061741_20276-no_umiCleanMerged.rds.gz",
            root / "GSM8061744_20276-no_CSP_counts_n1000.txt.gz",
        ),
        (
            "cytokine_treated",
            root / "GSM8061742_20276plus_umiCleanMerged.rds.gz",
            root / "GSM8061745_20276plus_CSP_counts_n1000.txt.gz",
        ),
    ]
    target_order = ["HLA_DR_COMPLEX", "PD_L1"]
    counts_parts = []
    expression_parts = []
    protein_parts = []
    cell_parts = []
    condition_parts = []
    for condition_index, (condition, rna_path, csp_path) in enumerate(configs):
        cache_path = cache_dir / f"gse254985_{condition}_rna_selected.tsv.gz"
        if not cache_path.exists():
            export_double_gzip_rds_selected(
                rna_path, gene_order=gene_order, output_path=cache_path
            )
        rna_cells, counts = read_feature_by_cell_tsv(
            cache_path, requested_features=gene_order
        )
        labels, protein_cells, protein_counts = read_all_feature_by_cell_tsv(csp_path)
        cells, counts, protein_counts = _align_cells(
            rna_cells, counts, protein_cells, protein_counts
        )
        if max_cells_per_arm > 0 and len(cells) > max_cells_per_arm:
            rng = np.random.default_rng(seed + condition_index)
            selected = np.sort(
                rng.choice(len(cells), size=max_cells_per_arm, replace=False)
            )
            cells = [cells[index] for index in selected]
            counts = counts[selected]
            protein_counts = protein_counts[selected]
        protein, targets = aggregate_target_adt(
            transform_target_adt_counts(protein_counts), labels, canonical_training_target
        )
        counts_parts.append(sparse.csr_matrix(counts))
        expression_parts.append(normalize_total_log1p(counts))
        protein_parts.append(_reorder_targets(protein, targets, target_order))
        cell_parts.append(np.asarray(cells, dtype=str))
        condition_parts.append(np.asarray([condition] * len(cells), dtype=str))
    return {
        "counts": sparse.vstack(counts_parts, format="csr"),
        "expression": np.vstack([_dense_float32(part) for part in expression_parts]),
        "protein": np.vstack(protein_parts),
        "targets": target_order,
        "cells": np.concatenate(cell_parts),
        "conditions": np.concatenate(condition_parts),
        "batches": np.concatenate(condition_parts),
    }


def fit_joint_scvi_latent(
    reference_counts,
    hgsoc_counts,
    *,
    reference_batches: Sequence[str],
    hgsoc_batches: Sequence[str],
    gene_order: Sequence[str],
    model_dir: Path,
    n_latent: int,
    max_epochs: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    import anndata as ad
    import scvi

    scvi.settings.seed = seed
    counts = sparse.vstack([reference_counts, hgsoc_counts], format="csr")
    adata = ad.AnnData(X=counts)
    adata.var_names = list(gene_order)
    adata.obs_names = [f"reference_{i}" for i in range(reference_counts.shape[0])] + [
        f"hgsoc_{i}" for i in range(hgsoc_counts.shape[0])
    ]
    adata.obs["batch"] = np.concatenate(
        [np.asarray(reference_batches, dtype=str), np.asarray(hgsoc_batches, dtype=str)]
    )
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    model = scvi.model.SCVI(
        adata,
        n_hidden=128,
        n_latent=n_latent,
        n_layers=1,
        dropout_rate=0.1,
        gene_likelihood="nb",
    )
    model.train(
        max_epochs=max_epochs,
        accelerator="cpu",
        devices=1,
        train_size=0.9,
        batch_size=256,
        early_stopping=True,
        enable_progress_bar=True,
        logger=False,
    )
    latent = model.get_latent_representation(batch_size=512)
    model_dir.parent.mkdir(parents=True, exist_ok=True)
    model.save(str(model_dir), overwrite=True, save_anndata=True)
    split = reference_counts.shape[0]
    return latent[:split].astype(np.float32), latent[split:].astype(np.float32)


def fit_reference_scvi_latent(
    reference_counts,
    *,
    gene_order: Sequence[str],
    model_dir: Path,
    n_latent: int,
    max_epochs: int,
    seed: int,
    progress: bool = True,
) -> np.ndarray:
    import anndata as ad
    import scvi

    if max_epochs < 1 or n_latent < 1:
        raise ValueError("Reference scVI epochs and latent dimension must be positive.")
    scvi.settings.seed = seed
    adata = ad.AnnData(X=sparse.csr_matrix(reference_counts))
    adata.var_names = list(gene_order)
    adata.obs_names = [f"reference_{index}" for index in range(adata.n_obs)]
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(
        adata,
        n_hidden=128,
        n_latent=n_latent,
        n_layers=1,
        dropout_rate=0.0,
        gene_likelihood="nb",
    )
    model.train(
        max_epochs=max_epochs,
        accelerator="cpu",
        devices=1,
        train_size=1.0,
        batch_size=min(256, adata.n_obs),
        early_stopping=False,
        enable_progress_bar=progress,
        logger=False,
    )
    latent = model.get_latent_representation(batch_size=512)
    model_dir.parent.mkdir(parents=True, exist_ok=True)
    model.save(str(model_dir), overwrite=True, save_anndata=True)
    return latent.astype(np.float32)


def validate_representation_manifest(manifest: Mapping[str, object]) -> None:
    required = {
        "mode",
        "fit_donors",
        "mapped_donors",
        "query_adaptation_epochs",
        "gene_order_sha256",
    }
    missing = required.difference(manifest)
    if missing:
        raise ValueError(f"Representation manifest is missing fields: {sorted(missing)}")
    fit_donors = {str(value) for value in manifest["fit_donors"]}
    mapped_donors = {str(value) for value in manifest["mapped_donors"]}
    if manifest["mode"] != "reference_only_inductive":
        raise ValueError("Strict benchmark representations must be reference-only inductive.")
    if not fit_donors or fit_donors.intersection(mapped_donors):
        raise ValueError("Representation fit and mapped donors must be non-overlapping.")
    if int(manifest["query_adaptation_epochs"]) != 0:
        raise ValueError("Strict benchmark queries cannot adapt the representation encoder.")
    if not str(manifest["gene_order_sha256"]):
        raise ValueError("Representation manifest requires a frozen gene-order hash.")


def map_scvi_query_latent(
    counts,
    *,
    batches: Sequence[str],
    gene_order: Sequence[str],
    reference_model_dir: Path,
    output_model_dir: Path,
    max_epochs: int,
    seed: int,
    progress: bool = True,
) -> np.ndarray:
    import anndata as ad
    import scvi

    scvi.settings.seed = seed
    query = ad.AnnData(X=sparse.csr_matrix(counts))
    query.var_names = list(gene_order)
    query.obs_names = [f"query_{i}" for i in range(query.n_obs)]
    query.obs["batch"] = np.asarray(batches, dtype=str)
    scvi.model.SCVI.prepare_query_anndata(query, str(reference_model_dir), inplace=True)
    model = scvi.model.SCVI.load_query_data(
        adata=query,
        reference_model=str(reference_model_dir),
        accelerator="cpu",
        device="auto",
    )
    if max_epochs < 0:
        raise ValueError("Query adaptation epochs cannot be negative.")
    if max_epochs > 0:
        model.train(
            max_epochs=max_epochs,
            accelerator="cpu",
            devices=1,
            train_size=1.0,
            batch_size=min(256, query.n_obs),
            early_stopping=False,
            plan_kwargs={"weight_decay": 0.0},
            enable_progress_bar=progress,
            logger=False,
        )
    else:
        # load_query_data restores reference weights but marks the wrapper untrained.
        model.is_trained = True
    latent = model.get_latent_representation(batch_size=512)
    model.save(str(output_model_dir), overwrite=True, save_anndata=True)
    return latent.astype(np.float32)


def _write_json(path: Path, value: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _prepared_cache_directory(args: argparse.Namespace) -> Path:
    seed = int(args.seeds[0])
    key = (
        f"g{args.max_genes}_r{args.max_ref_cells_per_donor}_"
        f"h{args.max_hgsoc_cells}_s{seed}"
    )
    return Path(args.output_dir) / "prepared" / key


def prepare_multidomain_data(args: argparse.Namespace) -> dict[str, object]:
    cache_dir = _prepared_cache_directory(args)
    manifest_path = cache_dir / "prepared_manifest.json"
    block_names = ("reference", "hgsoc", "gse128639", "gse254985")
    if manifest_path.exists() and all(
        (cache_dir / f"{name}_counts.npz").exists()
        and (cache_dir / f"{name}_arrays.npz").exists()
        for name in block_names
    ):
        with open(manifest_path, encoding="utf-8") as handle:
            manifest = json.load(handle)
        return {
            "cache_dir": cache_dir,
            "manifest": manifest,
            **{name: load_prepared_block(cache_dir, name) for name in block_names},
            "gene_embedding": np.load(cache_dir / "gene_embedding.npy", mmap_mode="r"),
            "protein_embedding": np.load(
                cache_dir / "protein_embedding.npy", mmap_mode="r"
            ),
        }

    data_root = Path(args.data_root)
    hgsoc_root = Path(args.hgsoc_root)
    reference_path = (
        data_root
        / "GSE194122"
        / "GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad"
    )
    gene_embedding_path = data_root / "embeddings" / "dna2vec_1w.npz"
    protein_embedding_path = data_root / "embeddings" / "per-protein.h5"
    identity_table = (
        Path(args.output_dir).parent
        / "tables"
        / "gse266577_eoc_identity_audit_by_cell.tsv.gz"
    )
    source_paths = [
        reference_path,
        gene_embedding_path,
        protein_embedding_path,
        hgsoc_root / "GSE266577_seurat_features.txt.gz",
        hgsoc_root / "GSE266577_barcodes.txt.gz",
        hgsoc_root / "GSE266577_metadata.txt.gz",
        hgsoc_root / "derived" / "eoc_csc" / "gse266577_eoc_counts_x.bin",
        hgsoc_root / "derived" / "eoc_csc" / "gse266577_eoc_counts_i.bin",
        hgsoc_root / "derived" / "eoc_csc" / "gse266577_eoc_counts_p.tsv",
        identity_table,
    ]
    missing_sources = [os.fspath(path) for path in source_paths if not path.exists()]
    if missing_sources:
        raise FileNotFoundError(f"Required multidomain sources are missing: {missing_sources}")

    reference_info = inspect_gse194122(reference_path)
    sampled_rows = sample_indices_by_group(
        reference_info["donors"],
        max_per_group=int(args.max_ref_cells_per_donor),
        seed=int(args.seeds[0]),
    )
    hgsoc_genes = set(
        read_line_names(hgsoc_root / "GSE266577_seurat_features.txt.gz")
    )
    gse128_genes = set(
        read_feature_names(data_root / "GSE128639" / "GSM3681518_MNC_RNA_counts.tsv.gz")
    )
    gse254_root = data_root / "GSE254985" / "extracted"
    gse254_genes = set(
        read_double_gzip_rds_gene_names(
            gse254_root / "GSM8061741_20276-no_umiCleanMerged.rds.gz"
        )
    ).intersection(
        read_double_gzip_rds_gene_names(
            gse254_root / "GSM8061742_20276plus_umiCleanMerged.rds.gz"
        )
    )
    embedding_genes = read_gene_embedding_names(gene_embedding_path)
    query_common = hgsoc_genes.intersection(gse128_genes, gse254_genes)
    reference_gene_names = list(reference_info["gene_names"])
    common_gene_order = [
        gene
        for gene in reference_gene_names
        if gene in query_common and gene in embedding_genes
    ]
    if len(common_gene_order) < int(args.max_genes):
        raise ValueError(
            f"Only {len(common_gene_order)} genes are shared across all benchmark arms."
        )

    reference_gene_index = {
        gene: int(index)
        for gene, index in zip(reference_gene_names, reference_info["gene_indices"])
    }
    target_adt = [
        (int(index), str(label), canonical_training_target(label))
        for index, label in zip(reference_info["adt_indices"], reference_info["adt_labels"])
        if canonical_training_target(label) in TARGET_ORDER
    ]
    target_columns = [item[0] for item in target_adt]
    combined = read_h5ad_csr_subset(
        reference_path,
        row_indices=sampled_rows,
        column_indices=[reference_gene_index[gene] for gene in common_gene_order]
        + target_columns,
        shape=tuple(reference_info["shape"]),
    )
    common_counts = combined[:, : len(common_gene_order)].tocsr()
    common_expression = normalize_total_log1p(common_counts)
    variances = sparse_feature_variance(common_expression)
    variance_scores = dict(zip(common_gene_order, variances))
    gene_order = select_model_genes(
        reference_genes=reference_gene_names,
        target_genes=query_common,
        embedding_genes=embedding_genes,
        variance_scores=variance_scores,
        forced_genes=FORCED_GENES,
        max_genes=int(args.max_genes),
    )
    common_positions = {gene: index for index, gene in enumerate(common_gene_order)}
    reference_counts = common_counts[:, [common_positions[gene] for gene in gene_order]]
    raw_target_counts = combined[:, len(common_gene_order) :].toarray()
    reference_protein, reference_targets = aggregate_target_adt(
        transform_target_adt_counts(raw_target_counts),
        [item[1] for item in target_adt],
        canonical_training_target,
    )
    reference_block = {
        "counts": reference_counts,
        "expression": _dense_float32(normalize_total_log1p(reference_counts)),
        "protein": _reorder_targets(
            reference_protein, reference_targets, TARGET_ORDER
        ),
        "targets": np.asarray(TARGET_ORDER),
        "cells": np.asarray(reference_info["cell_names"])[sampled_rows],
        "donors": np.asarray(reference_info["donors"])[sampled_rows],
        "batches": np.asarray(reference_info["batches"])[sampled_rows],
    }
    hgsoc_block = load_hgsoc_selected_counts(
        hgsoc_root,
        gene_order=gene_order,
        max_cells=int(args.max_hgsoc_cells),
        seed=int(args.seeds[0]),
        identity_table=identity_table,
    )
    gse128_block = load_gse128639_query(
        data_root,
        gene_order=gene_order,
        max_cells=max(500, int(args.max_ref_cells_per_donor)),
        seed=int(args.seeds[0]),
    )
    gse254_block = load_gse254985_query(
        data_root,
        gene_order=gene_order,
        cache_dir=cache_dir / "rds_exports",
        max_cells_per_arm=max(500, int(args.max_ref_cells_per_donor)),
        seed=int(args.seeds[0]),
    )
    gene_embedding = load_gene_embeddings(gene_embedding_path, gene_order)
    protein_embedding = load_protein_embeddings(
        protein_embedding_path,
        {target: TARGET_PROTEIN_ENTRIES[target] for target in TARGET_ORDER},
    )

    block_manifests = {
        "reference": save_prepared_block(cache_dir, "reference", reference_block),
        "hgsoc": save_prepared_block(cache_dir, "hgsoc", hgsoc_block),
        "gse128639": save_prepared_block(cache_dir, "gse128639", gse128_block),
        "gse254985": save_prepared_block(cache_dir, "gse254985", gse254_block),
    }
    np.save(cache_dir / "gene_embedding.npy", gene_embedding, allow_pickle=False)
    np.save(cache_dir / "protein_embedding.npy", protein_embedding, allow_pickle=False)
    gene_order_sha256 = hashlib.sha256("\n".join(gene_order).encode("utf-8")).hexdigest()
    manifest = {
        "schema_version": 1,
        "reference_dataset": "GSE194122_CITE",
        "reference_sampling": "donor_balanced",
        "reference_donor_count": int(np.unique(reference_block["donors"]).size),
        "gene_order": gene_order,
        "gene_order_sha256": gene_order_sha256,
        "target_order": list(TARGET_ORDER),
        "adt_transform": "targetwise_log1p_raw_counts",
        "feature_ranking_fit_scope": "reference_only",
        "blocks": block_manifests,
        "embedding_sha256": {
            "gene_embedding.npy": sha256_file(cache_dir / "gene_embedding.npy"),
            "protein_embedding.npy": sha256_file(cache_dir / "protein_embedding.npy"),
        },
        "source_sha256": {
            os.fspath(path): sha256_file(path) for path in source_paths
        },
        "parameters": {
            "max_genes": int(args.max_genes),
            "max_ref_cells_per_donor": int(args.max_ref_cells_per_donor),
            "max_hgsoc_cells": int(args.max_hgsoc_cells),
            "seed": int(args.seeds[0]),
        },
    }
    _write_json(manifest_path, manifest)
    return {
        "cache_dir": cache_dir,
        "manifest": manifest,
        "reference": reference_block,
        "hgsoc": hgsoc_block,
        "gse128639": gse128_block,
        "gse254985": gse254_block,
        "gene_embedding": gene_embedding,
        "protein_embedding": protein_embedding,
    }


def _stored_npy_layout(
    path: os.PathLike[str] | str,
    member_name: str,
) -> tuple[tuple[int, ...], bool, np.dtype, int]:
    from numpy.lib import format as npy_format

    path = os.fspath(path)
    with zipfile.ZipFile(path) as archive:
        info = archive.getinfo(member_name)
        if info.compress_type != zipfile.ZIP_STORED:
            raise ValueError(f"NPZ member must be uncompressed for mmap: {member_name}")
        header_offset = info.header_offset
    with open(path, "rb") as handle:
        handle.seek(header_offset)
        local_header = handle.read(30)
        if len(local_header) != 30 or local_header[:4] != b"PK\x03\x04":
            raise ValueError(f"Invalid ZIP header for {member_name}")
        filename_length, extra_length = struct.unpack("<HH", local_header[26:30])
        handle.seek(header_offset + 30 + filename_length + extra_length)
        version = npy_format.read_magic(handle)
        if version == (1, 0):
            shape, fortran_order, dtype = npy_format.read_array_header_1_0(handle)
        else:
            shape, fortran_order, dtype = npy_format.read_array_header_2_0(handle)
        data_offset = handle.tell()
    return tuple(shape), bool(fortran_order), np.dtype(dtype), data_offset


def load_gene_embeddings(
    path: os.PathLike[str] | str,
    gene_order: Sequence[str],
) -> np.ndarray:
    path = os.fspath(path)
    gene_shape, gene_fortran, gene_dtype, gene_offset = _stored_npy_layout(
        path, "gene.npy"
    )
    embedding_shape, embedding_fortran, embedding_dtype, embedding_offset = (
        _stored_npy_layout(path, "embedding.npy")
    )
    if gene_fortran or len(gene_shape) != 1 or embedding_shape[0] != gene_shape[0]:
        raise ValueError("Gene and embedding arrays are not aligned.")
    genes = np.memmap(
        path,
        mode="r",
        dtype=gene_dtype,
        offset=gene_offset,
        shape=gene_shape,
        order="C",
    )
    index = {str(gene): position for position, gene in enumerate(genes)}
    missing = [gene for gene in gene_order if gene not in index]
    if missing:
        raise KeyError(f"Missing official gene embeddings: {missing}")
    embeddings = np.memmap(
        path,
        mode="r",
        dtype=embedding_dtype,
        offset=embedding_offset,
        shape=embedding_shape,
        order="F" if embedding_fortran else "C",
    )
    positions = np.asarray([index[gene] for gene in gene_order], dtype=int)
    return np.asarray(embeddings[positions], dtype=np.float32)


def load_protein_embeddings(
    path: os.PathLike[str] | str,
    entries_by_target: Mapping[str, Sequence[str]],
) -> np.ndarray:
    import h5py

    rows = []
    with h5py.File(path, "r") as handle:
        for target, entries in entries_by_target.items():
            missing = [entry for entry in entries if entry not in handle]
            if missing:
                raise KeyError(f"Missing protein embeddings for {target}: {missing}")
            rows.append(
                np.mean(
                    np.stack(
                        [np.asarray(handle[entry], dtype=np.float32) for entry in entries],
                        axis=0,
                    ),
                    axis=0,
                    dtype=np.float32,
                )
            )
    if not rows:
        raise ValueError("At least one protein target is required.")
    return np.stack(rows, axis=0).astype(np.float32, copy=False)


def mean_chain_embedding(
    embeddings_by_entry: Mapping[str, np.ndarray],
    entries: Sequence[str],
) -> np.ndarray:
    if not entries:
        raise ValueError("At least one protein chain is required.")
    missing = [entry for entry in entries if entry not in embeddings_by_entry]
    if missing:
        raise KeyError(f"Missing protein embedding entries: {missing}")
    arrays = [np.asarray(embeddings_by_entry[entry], dtype=np.float32) for entry in entries]
    dimensions = {array.shape for array in arrays}
    if len(dimensions) != 1:
        raise ValueError("Protein-chain embeddings have inconsistent dimensions.")
    return np.mean(np.stack(arrays, axis=0), axis=0, dtype=np.float32)


def align_feature_matrix(
    matrix,
    source_names: Sequence[str],
    target_order: Sequence[str],
):
    source_names = list(source_names)
    if len(source_names) != len(set(source_names)):
        raise ValueError("Source feature names must be unique.")
    index = {name: position for position, name in enumerate(source_names)}
    missing = [name for name in target_order if name not in index]
    if missing:
        raise KeyError(f"Missing required features: {missing}")
    positions = [index[name] for name in target_order]
    if sparse.issparse(matrix):
        return matrix[:, positions]
    return np.asarray(matrix)[:, positions]


def factorized_sequence_cell_projection(
    layer: nn.Linear,
    sequence_embedding: torch.Tensor,
    cell_embedding: torch.Tensor,
) -> torch.Tensor:
    if sequence_embedding.ndim != 2 or cell_embedding.ndim != 2:
        raise ValueError("Factorized projection requires shared sequence and batched cell matrices.")
    sequence_dim = sequence_embedding.shape[1]
    if sequence_dim + cell_embedding.shape[1] != layer.in_features:
        raise ValueError("Sequence and cell dimensions do not align with the linear layer.")
    sequence_projection = F.linear(
        sequence_embedding,
        layer.weight[:, :sequence_dim],
        layer.bias,
    )
    cell_projection = F.linear(
        cell_embedding,
        layer.weight[:, sequence_dim:],
        bias=None,
    )
    return sequence_projection.unsqueeze(0) + cell_projection.unsqueeze(1)


class ProTrans(nn.Module):
    """Configurable reference-train/target-infer form of the official model."""

    def __init__(
        self,
        protein_dim: int = 1024,
        gene_dim: int = 10000,
        cell_dim: int = 100,
        hidden_dim: int = 500,
        heads: int = 10,
        head_dim: int = 50,
    ) -> None:
        super().__init__()
        if hidden_dim != heads * head_dim:
            raise ValueError("hidden_dim must equal heads * head_dim.")
        self.config = {
            "protein_dim": int(protein_dim),
            "gene_dim": int(gene_dim),
            "cell_dim": int(cell_dim),
            "hidden_dim": int(hidden_dim),
            "heads": int(heads),
            "head_dim": int(head_dim),
        }
        self.heads = heads
        self.head_dim = head_dim
        self.query_linear = nn.Linear(protein_dim + cell_dim, hidden_dim)
        self.key_linear = nn.Linear(gene_dim + cell_dim, hidden_dim)
        self.linear_1 = nn.Linear(heads, 10)
        self.linear_2 = nn.Linear(10, 1)

    def forward(
        self,
        protein_embedding: torch.Tensor,
        gene_embedding: torch.Tensor,
        cell_embedding: torch.Tensor,
        expression: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        shared_protein = protein_embedding.ndim == 2
        shared_gene = gene_embedding.ndim == 2
        if protein_embedding.ndim not in (2, 3) or gene_embedding.ndim not in (2, 3):
            raise ValueError("Sequence embeddings must be shared 2D or batched 3D tensors.")
        batch_size = expression.shape[0]
        protein_count = protein_embedding.shape[-2]
        gene_count = gene_embedding.shape[-2]
        if (
            (not shared_protein and protein_embedding.shape[0] != batch_size)
            or (not shared_gene and gene_embedding.shape[0] != batch_size)
            or expression.shape != (batch_size, gene_count)
        ):
            raise ValueError("Batch or gene dimensions are inconsistent.")
        if cell_embedding.shape[0] != batch_size:
            raise ValueError("Cell embedding batch dimension is inconsistent.")

        if shared_protein:
            query = factorized_sequence_cell_projection(
                self.query_linear, protein_embedding, cell_embedding
            )
        else:
            protein_cell = cell_embedding[:, None, :].expand(-1, protein_count, -1)
            query = self.query_linear(
                torch.cat([protein_embedding, protein_cell], dim=2)
            )
        if shared_gene:
            key = factorized_sequence_cell_projection(
                self.key_linear, gene_embedding, cell_embedding
            )
        else:
            gene_cell = cell_embedding[:, None, :].expand(-1, gene_count, -1)
            key = self.key_linear(torch.cat([gene_embedding, gene_cell], dim=2))
        query = query.view(batch_size, protein_count, self.heads, self.head_dim)
        key = key.view(batch_size, gene_count, self.heads, self.head_dim)
        query = query.permute(0, 2, 1, 3)
        key = key.permute(0, 2, 1, 3)

        attention = torch.einsum("bhpd,bhgd->bhpg", query, key)
        attention = torch.softmax(attention / math.sqrt(self.head_dim), dim=3)
        protein_values = torch.einsum("bhpg,bg->bhp", attention, expression)
        protein_values = protein_values.permute(0, 2, 1)
        prediction = self.linear_2(F.relu(self.linear_1(protein_values)))
        prediction = F.relu(prediction).squeeze(2)
        return prediction, attention.sum(dim=1)


def save_model_bundle(
    path: os.PathLike[str] | str,
    model: ProTrans,
    metadata: Mapping[str, object],
) -> None:
    required = {"gene_order", "protein_order", "normalization", "seed", "model_config"}
    missing = required.difference(metadata)
    if missing:
        raise ValueError(f"Model metadata is missing required fields: {sorted(missing)}")
    torch.save(
        {"state_dict": model.state_dict(), "metadata": dict(metadata)},
        os.fspath(path),
    )


def load_model_bundle(
    path: os.PathLike[str] | str,
) -> tuple[ProTrans, dict[str, object]]:
    bundle = torch.load(os.fspath(path), map_location="cpu", weights_only=False)
    metadata = dict(bundle["metadata"])
    model = ProTrans(**metadata["model_config"])
    model.load_state_dict(bundle["state_dict"], strict=True)
    model.eval()
    return model, metadata


def _dense_float32(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float32)


def _set_random_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def initialize_output_bias(model: ProTrans, protein_truth) -> None:
    values = np.asarray(protein_truth, dtype=float)
    if values.size == 0 or not np.isfinite(values).all():
        raise ValueError("Output-bias initialization requires finite protein truth.")
    initial_bias = max(0.1, float(np.median(values)))
    with torch.no_grad():
        model.linear_2.bias.fill_(initial_bias)


def train_protrans(
    expression,
    protein_truth,
    cell_latent,
    gene_embedding,
    protein_embedding,
    *,
    model_config: Mapping[str, int],
    seed: int,
    epochs: int = 50,
    batch_size: int = 32,
    learning_rate: float = 1e-3,
    patience: int = 10,
    device: str = "cpu",
) -> tuple[ProTrans, list[float]]:
    expression = _dense_float32(expression)
    protein_truth = _dense_float32(protein_truth)
    cell_latent = _dense_float32(cell_latent)
    gene_embedding = _dense_float32(gene_embedding)
    protein_embedding = _dense_float32(protein_embedding)
    if expression.shape[0] != protein_truth.shape[0] or expression.shape[0] != cell_latent.shape[0]:
        raise ValueError("Expression, protein truth, and cell latent must share cells.")
    if expression.shape[1] != gene_embedding.shape[0]:
        raise ValueError("Expression genes must align with gene embeddings.")
    if protein_truth.shape[1] != protein_embedding.shape[0]:
        raise ValueError("Protein truth must align with protein embeddings.")
    if epochs < 1 or batch_size < 1 or patience < 1:
        raise ValueError("epochs, batch_size, and patience must be positive.")

    _set_random_seed(seed)
    model = ProTrans(**dict(model_config)).to(device)
    initialize_output_bias(model, protein_truth)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    criterion = nn.MSELoss()
    dataset = torch.utils.data.TensorDataset(
        torch.from_numpy(expression),
        torch.from_numpy(protein_truth),
        torch.from_numpy(cell_latent),
    )
    generator = torch.Generator().manual_seed(seed)
    loader = torch.utils.data.DataLoader(
        dataset,
        batch_size=min(batch_size, len(dataset)),
        shuffle=True,
        generator=generator,
        drop_last=False,
    )
    gene_tensor = torch.from_numpy(gene_embedding).to(device)
    protein_tensor = torch.from_numpy(protein_embedding).to(device)
    history: list[float] = []
    best_loss = float("inf")
    best_state = None
    stale_epochs = 0
    for _ in range(epochs):
        model.train()
        total_loss = 0.0
        total_cells = 0
        for batch_expression, batch_protein, batch_latent in loader:
            batch_expression = batch_expression.to(device)
            batch_protein = batch_protein.to(device)
            batch_latent = batch_latent.to(device)
            optimizer.zero_grad(set_to_none=True)
            prediction, _ = model(
                protein_tensor,
                gene_tensor,
                batch_latent,
                batch_expression,
            )
            loss = criterion(prediction, batch_protein)
            loss.backward()
            optimizer.step()
            total_loss += float(loss.detach().cpu()) * batch_expression.shape[0]
            total_cells += batch_expression.shape[0]
        epoch_loss = total_loss / total_cells
        history.append(epoch_loss)
        if epoch_loss < best_loss - 1e-8:
            best_loss = epoch_loss
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }
            stale_epochs = 0
        else:
            stale_epochs += 1
            if stale_epochs >= patience:
                break
    if best_state is None:
        raise RuntimeError("Training did not produce a finite model state.")
    model = model.to("cpu")
    model.load_state_dict(best_state, strict=True)
    model.eval()
    return model, history


def predict_protrans(
    model: ProTrans,
    expression,
    cell_latent,
    gene_embedding,
    protein_embedding,
    *,
    batch_size: int = 64,
    device: str = "cpu",
) -> np.ndarray:
    expression = _dense_float32(expression)
    cell_latent = _dense_float32(cell_latent)
    gene_embedding = _dense_float32(gene_embedding)
    protein_embedding = _dense_float32(protein_embedding)
    if expression.shape[0] != cell_latent.shape[0]:
        raise ValueError("Expression and cell latent must share cells.")
    model = model.to(device)
    model.eval()
    gene_tensor = torch.from_numpy(gene_embedding).to(device)
    protein_tensor = torch.from_numpy(protein_embedding).to(device)
    predictions = []
    with torch.no_grad():
        for start in range(0, expression.shape[0], batch_size):
            stop = min(start + batch_size, expression.shape[0])
            batch_expression = torch.from_numpy(expression[start:stop]).to(device)
            batch_latent = torch.from_numpy(cell_latent[start:stop]).to(device)
            prediction, _ = model(
                protein_tensor,
                gene_tensor,
                batch_latent,
                batch_expression,
            )
            predictions.append(prediction.cpu().numpy())
    model.to("cpu")
    return np.concatenate(predictions, axis=0).astype(np.float32, copy=False)


def fit_ridge_baseline(expression, protein_truth, *, alpha: float = 1.0) -> Ridge:
    expression = _dense_float32(expression)
    protein_truth = _dense_float32(protein_truth)
    if expression.shape[0] != protein_truth.shape[0]:
        raise ValueError("Expression and protein truth must share cells.")
    return Ridge(alpha=alpha).fit(expression, protein_truth)


def predict_ridge_baseline(model: Ridge, expression) -> np.ndarray:
    prediction = np.asarray(model.predict(_dense_float32(expression)), dtype=np.float32)
    if prediction.ndim == 1:
        prediction = prediction[:, None]
    return prediction


def cognate_rna_score(
    expression,
    *,
    gene_order: Sequence[str],
    target: str,
) -> np.ndarray:
    if target not in TARGET_COGNATE_GENES:
        raise KeyError(f"No cognate RNA definition is registered for {target}.")
    index = {gene: position for position, gene in enumerate(gene_order)}
    positions = [index[gene] for gene in TARGET_COGNATE_GENES[target] if gene in index]
    if not positions:
        raise KeyError(f"No cognate genes for {target} occur in the frozen gene order.")
    values = expression[:, positions]
    if sparse.issparse(values):
        return np.asarray(values.mean(axis=1)).reshape(-1)
    return np.mean(np.asarray(values, dtype=float), axis=1)


def _safe_correlation(function, truth: np.ndarray, prediction: np.ndarray) -> float:
    if np.std(truth) == 0 or np.std(prediction) == 0:
        return float("nan")
    return float(function(truth, prediction).statistic)


def compute_prediction_metrics(
    truth: Sequence[float],
    prediction: Sequence[float],
) -> dict[str, float]:
    truth = np.asarray(truth, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    if truth.shape != prediction.shape or truth.ndim != 1 or truth.size < 2:
        raise ValueError("Truth and prediction must be aligned one-dimensional vectors.")
    if not np.isfinite(truth).all() or not np.isfinite(prediction).all():
        raise ValueError("Metrics require finite values.")
    cosine_denominator = float(np.linalg.norm(truth) * np.linalg.norm(prediction))
    if np.std(prediction) == 0:
        intercept = float("nan")
        slope = float("nan")
    else:
        design = np.column_stack([np.ones(prediction.size), prediction])
        intercept, slope = np.linalg.lstsq(design, truth, rcond=None)[0]
    cosine = (
        float(np.dot(truth, prediction) / cosine_denominator)
        if cosine_denominator
        else float("nan")
    )
    return {
        "spearman": _safe_correlation(spearmanr, truth, prediction),
        "pearson": _safe_correlation(pearsonr, truth, prediction),
        "mae": float(mean_absolute_error(truth, prediction)),
        "cosine": cosine,
        "calibration_intercept": float(intercept),
        "calibration_slope": float(slope),
    }


def compute_group_metric_records(
    truth: Sequence[float],
    prediction: Sequence[float],
    groups: Sequence[str],
    *,
    metric: str = "spearman",
) -> list[dict[str, object]]:
    truth = np.asarray(truth, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    groups = np.asarray(groups, dtype=str)
    if truth.shape != prediction.shape or truth.shape != groups.shape or truth.ndim != 1:
        raise ValueError("Truth, prediction, and groups must be aligned vectors.")
    records = []
    for group in np.unique(groups):
        group_mask = groups == group
        metrics = compute_prediction_metrics(truth[group_mask], prediction[group_mask])
        if metric not in metrics:
            raise ValueError(f"Unsupported metric: {metric}")
        records.append(
            {
                "group": str(group),
                "cell_count": int(np.sum(group_mask)),
                "value": float(metrics[metric]),
            }
        )
    return records


def cluster_bootstrap_metric(
    truth: Sequence[float],
    prediction: Sequence[float],
    groups: Sequence[str],
    *,
    metric: str = "spearman",
    iterations: int = 2000,
    seed: int = 260716,
) -> dict[str, float | int]:
    truth = np.asarray(truth, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    groups = np.asarray(groups, dtype=str)
    if truth.shape != prediction.shape or truth.shape != groups.shape or truth.ndim != 1:
        raise ValueError("Truth, prediction, and groups must be aligned vectors.")
    if iterations < 1 or np.unique(groups).size < 3:
        raise ValueError("Cluster bootstrap requires at least three groups and one iteration.")
    if metric not in compute_prediction_metrics(truth, prediction):
        raise ValueError(f"Unsupported metric: {metric}")
    unique_groups = np.unique(groups)
    group_metrics = []
    for group in unique_groups:
        group_mask = groups == group
        value = compute_prediction_metrics(truth[group_mask], prediction[group_mask])[metric]
        if np.isfinite(value):
            group_metrics.append(float(value))
    group_metrics = np.asarray(group_metrics, dtype=float)
    if group_metrics.size < 3:
        raise ValueError("Fewer than three groups have a finite metric.")
    rng = np.random.default_rng(seed)
    draws = np.mean(
        rng.choice(group_metrics, size=(iterations, group_metrics.size), replace=True),
        axis=1,
    )
    return {
        "estimate": float(np.mean(group_metrics)),
        "ci_low": float(np.quantile(draws, 0.025)),
        "ci_high": float(np.quantile(draws, 0.975)),
        "iterations": int(iterations),
        "finite_draws": int(np.isfinite(draws).sum()),
        "donor_count": int(unique_groups.size),
        "finite_donor_count": int(group_metrics.size),
    }


def fit_knn_ood(
    reference_latent: np.ndarray,
    reference_groups: Sequence[str] | None = None,
    k: int = 10,
    quantile: float = 0.95,
) -> dict[str, object]:
    reference = np.asarray(reference_latent, dtype=np.float32)
    if reference.ndim != 2 or reference.shape[0] <= k or k < 1:
        raise ValueError("Reference latent matrix must contain more than k cells.")
    if not 0 < quantile < 1:
        raise ValueError("OOD quantile must lie strictly between zero and one.")
    if reference_groups is None:
        neighbors = NearestNeighbors(n_neighbors=k + 1).fit(reference)
        distances, _ = neighbors.kneighbors(reference)
        calibration_distances = distances[:, k]
        calibration = "leave_one_cell_out"
    else:
        groups = np.asarray(reference_groups, dtype=str)
        if groups.shape != (reference.shape[0],) or np.unique(groups).size < 2:
            raise ValueError("Reference groups must align with cells and contain at least two groups.")
        calibration_distances = []
        for index, group in enumerate(groups):
            candidates = reference[groups != group]
            if candidates.shape[0] < k:
                raise ValueError("Each group holdout must leave at least k reference cells.")
            neighbors = NearestNeighbors(n_neighbors=k).fit(candidates)
            distances, _ = neighbors.kneighbors(reference[index : index + 1])
            calibration_distances.append(float(distances[0, k - 1]))
        calibration_distances = np.asarray(calibration_distances, dtype=float)
        calibration = "leave_one_group_out"
    threshold = float(np.quantile(calibration_distances, quantile))
    return {
        "reference_latent": reference,
        "k": int(k),
        "quantile": float(quantile),
        "threshold": threshold,
        "calibration": calibration,
    }


def score_knn_ood(
    model: Mapping[str, object],
    target_latent: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    reference = np.asarray(model["reference_latent"], dtype=np.float32)
    target = np.asarray(target_latent, dtype=np.float32)
    k = int(model["k"])
    neighbors = NearestNeighbors(n_neighbors=k).fit(reference)
    distances, _ = neighbors.kneighbors(target)
    scores = distances[:, k - 1]
    return scores, scores <= float(model["threshold"])


def patient_equal_in_domain_fraction(
    in_domain: Sequence[bool],
    patients: Sequence[str],
) -> dict[str, object]:
    flags = np.asarray(in_domain, dtype=bool)
    patients = np.asarray(patients, dtype=str)
    if flags.shape != patients.shape or flags.ndim != 1 or flags.size == 0:
        raise ValueError("In-domain flags and patient labels must be aligned vectors.")
    patient_fractions = {
        patient: float(np.mean(flags[patients == patient])) for patient in np.unique(patients)
    }
    return {
        "cell_fraction": float(np.mean(flags)),
        "patient_equal_fraction": float(np.mean(list(patient_fractions.values()))),
        "patient_count": len(patient_fractions),
        "patient_fractions": patient_fractions,
    }


def compute_patient_paired_effect(
    prediction: Sequence[float],
    patients: Sequence[str],
    stages: Sequence[str],
    *,
    pre_label: str,
    post_label: str,
) -> dict[str, object]:
    prediction = np.asarray(prediction, dtype=float)
    patients = np.asarray(patients, dtype=str)
    stages = np.asarray(stages, dtype=str)
    if prediction.shape != patients.shape or prediction.shape != stages.shape:
        raise ValueError("Predictions, patients, and stages must be aligned vectors.")
    paired_patients = []
    effects = []
    for patient in np.unique(patients):
        patient_mask = patients == patient
        pre = prediction[patient_mask & (stages == pre_label)]
        post = prediction[patient_mask & (stages == post_label)]
        if pre.size and post.size:
            paired_patients.append(patient)
            effects.append(float(np.mean(post) - np.mean(pre)))
    if not effects:
        raise ValueError("No patients contain both requested treatment stages.")
    effect_array = np.asarray(effects, dtype=float)
    return {
        "patient_count": len(paired_patients),
        "patients": paired_patients,
        "patient_effects": effect_array,
        "mean_effect": float(np.mean(effect_array)),
        "median_effect": float(np.median(effect_array)),
        "patient_effect_sd": float(np.std(effect_array, ddof=1)) if effect_array.size > 1 else 0.0,
    }


def read_eligible_patients(
    path: os.PathLike[str] | str,
    column: str = "eoc_pair_min10",
) -> set[str]:
    eligible = set()
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get(column, "").upper() == "TRUE":
                eligible.add(row["patient_id"])
    if not eligible:
        raise ValueError(f"No eligible patients were found in {path} using {column}.")
    return eligible


def not_worse_than_baseline(
    *,
    model_metric: float,
    baseline_metric: float,
    margin: float = 0.05,
) -> bool:
    if margin < 0:
        raise ValueError("The non-inferiority margin must be non-negative.")
    if not np.isfinite(model_metric) or not np.isfinite(baseline_metric):
        return False
    return bool(model_metric >= baseline_metric - margin)


def paired_noninferiority_summary(
    model_metrics: Sequence[float],
    baseline_metrics: Sequence[float],
    *,
    margin: float = 0.05,
    iterations: int = 2000,
    seed: int = 260716,
) -> dict[str, object]:
    model = np.asarray(model_metrics, dtype=float)
    baseline = np.asarray(baseline_metrics, dtype=float)
    if model.shape != baseline.shape or model.ndim != 1 or model.size < 3:
        raise ValueError("Paired non-inferiority requires at least three aligned units.")
    if margin < 0 or iterations < 1:
        raise ValueError("Margin must be non-negative and iterations must be positive.")
    if not np.isfinite(model).all() or not np.isfinite(baseline).all():
        raise ValueError("Paired non-inferiority metrics must be finite.")
    differences = model - baseline
    rng = np.random.default_rng(seed)
    draws = np.mean(
        rng.choice(differences, size=(iterations, differences.size), replace=True),
        axis=1,
    )
    ci_low = float(np.quantile(draws, 0.025))
    return {
        "unit_count": int(differences.size),
        "mean_difference": float(np.mean(differences)),
        "ci_low": ci_low,
        "ci_high": float(np.quantile(draws, 0.975)),
        "margin": float(margin),
        "iterations": int(iterations),
        "passed": bool(ci_low > -margin),
    }


def summarize_seed_patient_effects(seed_patient_effects) -> dict[str, object]:
    effects = np.asarray(seed_patient_effects, dtype=float)
    if effects.ndim != 2 or effects.shape[0] != 5 or effects.shape[1] < 2:
        raise ValueError("Seed stability requires a 5-by-at-least-2 patient matrix.")
    if not np.isfinite(effects).all():
        raise ValueError("Seed-by-patient effects must be finite.")
    seed_means = np.mean(effects, axis=1)
    patient_means = np.mean(effects, axis=0)
    return {
        "seed_count": int(effects.shape[0]),
        "patient_count": int(effects.shape[1]),
        "seed_means": seed_means.tolist(),
        "patient_means": patient_means.tolist(),
        "all_seed_means_positive": bool(np.all(seed_means > 0)),
        "seed_effect_sd": float(np.std(seed_means, ddof=1)),
        "patient_effect_sd": float(np.std(patient_means, ddof=1)),
    }


def evaluate_gate3_critical_stage(
    *,
    donor_summary: Mapping[str, object],
    independent_model_metrics: Sequence[float],
    independent_baseline_metrics: Sequence[float],
    hgsoc_patient_equal_in_domain_fraction: float,
    bundle_audit: Mapping[str, object],
    noninferiority_margin: float = 0.05,
    bootstrap_iterations: int = 2000,
    seed: int = 260716,
) -> dict[str, object]:
    required_donor = {
        "estimate",
        "ci_low",
        "donor_count",
        "finite_donor_count",
        "finite_draws",
        "iterations",
    }
    missing = required_donor.difference(donor_summary)
    if missing:
        raise ValueError(f"Donor summary is missing fields: {sorted(missing)}")
    donor_count = int(donor_summary["donor_count"])
    finite_donor_count = int(donor_summary["finite_donor_count"])
    finite_draws = int(donor_summary["finite_draws"])
    iterations = int(donor_summary["iterations"])
    if (
        donor_count < 1
        or finite_donor_count < 0
        or finite_donor_count > donor_count
        or iterations < 1
        or finite_draws < 0
        or finite_draws > iterations
        or not 0 <= hgsoc_patient_equal_in_domain_fraction <= 1
    ):
        raise ValueError("Critical-stage Gate 3 inputs contain invalid ranges.")
    independent_model = np.asarray(independent_model_metrics, dtype=float)
    noninferiority = paired_noninferiority_summary(
        independent_model,
        independent_baseline_metrics,
        margin=noninferiority_margin,
        iterations=bootstrap_iterations,
        seed=seed,
    )
    required_bundle = {
        "roundtrip_match",
        "artifact_hashes_match",
        "encoder_present",
        "split_manifest_present",
    }
    conditions = {
        "donor_holdout_performance": bool(
            float(donor_summary["estimate"]) >= 0.30
            and float(donor_summary["ci_low"]) > 0
            and donor_count >= 5
            and finite_donor_count / donor_count >= 0.80
            and finite_draws / iterations >= 0.80
        ),
        "independent_transfer": bool(
            np.isfinite(independent_model).all()
            and np.all(independent_model > 0)
            and noninferiority["passed"]
        ),
        "hgsoc_ood_coverage": bool(
            hgsoc_patient_equal_in_domain_fraction >= 0.70
        ),
        "reproducible_bundle": bool(
            required_bundle.issubset(bundle_audit)
            and all(bundle_audit.get(key) is True for key in required_bundle)
        ),
    }
    failed = [name for name, passed in conditions.items() if not passed]
    return {
        "proceed_to_five_seed": not failed,
        "failed_conditions": failed,
        "independent_noninferiority": noninferiority,
        **conditions,
    }


def evaluate_gate3(
    *,
    donor_summary: Mapping[str, object],
    independent_model_metrics: Sequence[float],
    independent_baseline_metrics: Sequence[float],
    hgsoc_in_domain: Sequence[bool],
    hgsoc_patients: Sequence[str],
    seed_patient_effects,
    bundle_audit: Mapping[str, object],
    noninferiority_margin: float = 0.05,
    bootstrap_iterations: int = 2000,
    seed: int = 260716,
) -> dict[str, object]:
    required_donor = {
        "estimate",
        "ci_low",
        "donor_count",
        "finite_donor_count",
        "finite_draws",
        "iterations",
    }
    missing = required_donor.difference(donor_summary)
    if missing:
        raise ValueError(f"Donor summary is missing fields: {sorted(missing)}")
    donor_spearman = float(donor_summary["estimate"])
    donor_ci_low = float(donor_summary["ci_low"])
    donor_count = int(donor_summary["donor_count"])
    finite_donor_count = int(donor_summary["finite_donor_count"])
    finite_draws = int(donor_summary["finite_draws"])
    iterations = int(donor_summary["iterations"])
    if (
        not np.isfinite(donor_spearman)
        or not np.isfinite(donor_ci_low)
        or donor_count < 1
        or finite_donor_count < 0
        or finite_donor_count > donor_count
        or iterations < 1
        or finite_draws < 0
        or finite_draws > iterations
    ):
        raise ValueError("Donor summary contains invalid values.")
    independent_model = np.asarray(independent_model_metrics, dtype=float)
    noninferiority = paired_noninferiority_summary(
        independent_model,
        independent_baseline_metrics,
        margin=noninferiority_margin,
        iterations=bootstrap_iterations,
        seed=seed,
    )
    ood = patient_equal_in_domain_fraction(hgsoc_in_domain, hgsoc_patients)
    stability = summarize_seed_patient_effects(seed_patient_effects)
    required_bundle = {
        "roundtrip_match",
        "artifact_hashes_match",
        "encoder_present",
        "split_manifest_present",
    }
    bundle_verified = required_bundle.issubset(bundle_audit) and all(
        bundle_audit.get(key) is True for key in required_bundle
    )
    conditions = {
        "donor_holdout_performance": bool(
            donor_spearman >= 0.30
            and donor_ci_low > 0
            and donor_count >= 5
            and finite_donor_count / donor_count >= 0.80
            and finite_draws / iterations >= 0.80
        ),
        "independent_transfer": bool(
            np.isfinite(independent_model).all()
            and np.all(independent_model > 0)
            and noninferiority["passed"]
        ),
        "hgsoc_ood_coverage": ood["patient_equal_fraction"] >= 0.70,
        "seed_direction_stability": stability["all_seed_means_positive"],
        "seed_variance_below_patient_variance": bool(
            stability["patient_effect_sd"] > 0
            and stability["seed_effect_sd"] < stability["patient_effect_sd"]
        ),
        "reproducible_bundle": bool(bundle_verified),
    }
    failed = [name for name, passed in conditions.items() if not passed]
    return {
        "gate3_passed": not failed,
        "failed_conditions": failed,
        "donor_summary": dict(donor_summary),
        "independent_noninferiority": noninferiority,
        "hgsoc_ood_summary": ood,
        "seed_stability_summary": stability,
        "bundle_audit": dict(bundle_audit),
        **conditions,
    }


def resolve_torch_device(requested: str) -> str:
    if requested != "auto":
        return requested
    if torch.cuda.is_available():
        return "cuda"
    if torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def _write_tsv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if not rows:
        raise ValueError(f"Cannot write an empty TSV: {path}")
    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def _benchmark_model_config(
    prepared: Mapping[str, object],
    *,
    hidden_dim: int,
    heads: int,
) -> dict[str, int]:
    if hidden_dim < 1 or heads < 1 or hidden_dim % heads:
        raise ValueError("Model hidden dimension must be positive and divisible by heads.")
    return {
        "protein_dim": int(np.asarray(prepared["protein_embedding"]).shape[1]),
        "gene_dim": int(np.asarray(prepared["gene_embedding"]).shape[1]),
        "cell_dim": int(prepared["scvi_latent"]),
        "hidden_dim": int(hidden_dim),
        "heads": int(heads),
        "head_dim": hidden_dim // heads,
    }


def _prediction_records(
    *,
    cells: Sequence[str],
    donors: Sequence[str],
    truth: np.ndarray,
    predictions: Mapping[str, np.ndarray],
    targets: Sequence[str],
    folds: Sequence[int] | None = None,
) -> list[dict[str, object]]:
    records = []
    for cell_index, (cell, donor) in enumerate(zip(cells, donors)):
        for target_index, target in enumerate(targets):
            record = {
                "cell": str(cell),
                "donor": str(donor),
                "target": str(target),
                "truth": float(truth[cell_index, target_index]),
            }
            if folds is not None:
                record["fold"] = int(folds[cell_index])
            for method, values in predictions.items():
                value = values[cell_index, target_index]
                record[method] = float(value) if np.isfinite(value) else ""
            records.append(record)
    return records


def run_donor_holdout_benchmark(
    prepared: Mapping[str, object],
    args: argparse.Namespace,
    *,
    seed: int,
    run_dir: Path,
) -> dict[str, object]:
    reference = prepared["reference"]
    donors = np.asarray(reference["donors"], dtype=str)
    targets = [str(value) for value in reference["targets"]]
    gene_order = list(prepared["manifest"]["gene_order"])
    gene_embedding = np.array(prepared["gene_embedding"], dtype=np.float32, copy=True)
    protein_embedding = np.array(
        prepared["protein_embedding"], dtype=np.float32, copy=True
    )
    model_config = _benchmark_model_config(
        {**prepared, "scvi_latent": int(args.scvi_latent)},
        hidden_dim=int(args.model_hidden_dim),
        heads=int(args.model_heads),
    )
    folds = make_group_folds(
        donors,
        n_splits=min(5, np.unique(donors).size),
        seed=seed,
    )
    prediction_by_method = {
        "scprotrans": np.full_like(reference["protein"], np.nan, dtype=np.float32),
        "ridge": np.full_like(reference["protein"], np.nan, dtype=np.float32),
        "cognate_rna": np.full_like(reference["protein"], np.nan, dtype=np.float32),
    }
    fold_assignment = np.full(len(donors), -1, dtype=int)
    representation_manifests = []
    device = resolve_torch_device(args.device)
    for fold_index, (train_index, test_index) in enumerate(folds, start=1):
        fold_dir = run_dir / "folds" / f"fold_{fold_index}"
        reference_model_dir = fold_dir / "scvi_reference"
        train_latent = fit_reference_scvi_latent(
            reference["counts"][train_index],
            gene_order=gene_order,
            model_dir=reference_model_dir,
            n_latent=int(args.scvi_latent),
            max_epochs=int(args.scvi_epochs),
            seed=seed + fold_index,
            progress=False,
        )
        test_latent = map_scvi_query_latent(
            reference["counts"][test_index],
            batches=["heldout"] * len(test_index),
            gene_order=gene_order,
            reference_model_dir=reference_model_dir,
            output_model_dir=fold_dir / "scvi_heldout",
            max_epochs=0,
            seed=seed + fold_index,
            progress=False,
        )
        manifest = {
            "mode": "reference_only_inductive",
            "fit_donors": sorted(np.unique(donors[train_index]).tolist()),
            "mapped_donors": sorted(np.unique(donors[test_index]).tolist()),
            "query_adaptation_epochs": 0,
            "gene_order_sha256": prepared["manifest"]["gene_order_sha256"],
        }
        validate_representation_manifest(manifest)
        _write_json(fold_dir / "representation_manifest.json", manifest)
        representation_manifests.append(manifest)

        model, _ = train_protrans(
            reference["expression"][train_index],
            reference["protein"][train_index],
            train_latent,
            gene_embedding,
            protein_embedding,
            model_config=model_config,
            seed=seed + fold_index,
            epochs=int(args.epochs),
            batch_size=int(args.batch_size),
            patience=int(args.patience),
            device=device,
        )
        prediction_by_method["scprotrans"][test_index] = predict_protrans(
            model,
            reference["expression"][test_index],
            test_latent,
            gene_embedding,
            protein_embedding,
            batch_size=int(args.batch_size),
            device=device,
        )
        ridge = fit_ridge_baseline(
            reference["expression"][train_index],
            reference["protein"][train_index],
            alpha=1.0,
        )
        prediction_by_method["ridge"][test_index] = predict_ridge_baseline(
            ridge, reference["expression"][test_index]
        )
        for target_index, target in enumerate(targets):
            try:
                prediction_by_method["cognate_rna"][test_index, target_index] = (
                    cognate_rna_score(
                        reference["expression"][test_index],
                        gene_order=gene_order,
                        target=target,
                    )
                )
            except KeyError:
                continue
        fold_assignment[test_index] = fold_index

    if fold_assignment.min() < 1 or not np.isfinite(
        prediction_by_method["scprotrans"]
    ).all():
        raise RuntimeError("Donor holdout did not produce complete out-of-fold predictions.")
    prediction_records = _prediction_records(
        cells=reference["cells"],
        donors=donors,
        truth=np.asarray(reference["protein"]),
        predictions=prediction_by_method,
        targets=targets,
        folds=fold_assignment,
    )
    _write_tsv(run_dir / "donor_holdout_predictions.tsv.gz", prediction_records)

    metric_rows = []
    for target_index, target in enumerate(targets):
        for method, values in prediction_by_method.items():
            if not np.isfinite(values[:, target_index]).all():
                continue
            for donor_record in compute_group_metric_records(
                reference["protein"][:, target_index],
                values[:, target_index],
                donors,
                metric="spearman",
            ):
                metric_rows.append(
                    {
                        "evaluation": "donor_holdout",
                        "target": target,
                        "method": method,
                        "unit": donor_record["group"],
                        "cell_count": donor_record["cell_count"],
                        "spearman": donor_record["value"],
                    }
                )
    _write_tsv(run_dir / "donor_holdout_metrics.tsv", metric_rows)
    hla_index = targets.index("HLA_DR_COMPLEX")
    hla_summary = cluster_bootstrap_metric(
        reference["protein"][:, hla_index],
        prediction_by_method["scprotrans"][:, hla_index],
        donors,
        metric="spearman",
        iterations=1000,
        seed=seed,
    )
    return {
        "hla_dr_summary": hla_summary,
        "metric_rows": metric_rows,
        "prediction_by_method": prediction_by_method,
        "fold_assignment": fold_assignment,
        "representation_manifests": representation_manifests,
        "model_config": model_config,
    }


def run_application_benchmark(
    prepared: Mapping[str, object],
    args: argparse.Namespace,
    *,
    seed: int,
    run_dir: Path,
) -> dict[str, object]:
    application_dir = run_dir / "application"
    application_dir.mkdir(parents=True, exist_ok=True)
    reference = prepared["reference"]
    gene_order = list(prepared["manifest"]["gene_order"])
    targets = [str(value) for value in reference["targets"]]
    gene_embedding = np.array(prepared["gene_embedding"], dtype=np.float32, copy=True)
    protein_embedding = np.array(
        prepared["protein_embedding"], dtype=np.float32, copy=True
    )
    model_config = _benchmark_model_config(
        {**prepared, "scvi_latent": int(args.scvi_latent)},
        hidden_dim=int(args.model_hidden_dim),
        heads=int(args.model_heads),
    )
    reference_scvi_dir = application_dir / "scvi_reference"
    reference_latent = fit_reference_scvi_latent(
        reference["counts"],
        gene_order=gene_order,
        model_dir=reference_scvi_dir,
        n_latent=int(args.scvi_latent),
        max_epochs=int(args.scvi_epochs),
        seed=seed,
        progress=False,
    )
    device = resolve_torch_device(args.device)
    model, history = train_protrans(
        reference["expression"],
        reference["protein"],
        reference_latent,
        gene_embedding,
        protein_embedding,
        model_config=model_config,
        seed=seed,
        epochs=int(args.epochs),
        batch_size=int(args.batch_size),
        patience=int(args.patience),
        device=device,
    )
    ridge = fit_ridge_baseline(
        reference["expression"], reference["protein"], alpha=1.0
    )
    representation_manifest = {
        "mode": "reference_only_inductive",
        "fit_donors": sorted(np.unique(reference["donors"]).tolist()),
        "mapped_donors": ["GSE128639", "GSE254985", "HGSOC_GSE266577"],
        "query_adaptation_epochs": 0,
        "gene_order_sha256": prepared["manifest"]["gene_order_sha256"],
    }
    validate_representation_manifest(representation_manifest)
    representation_path = application_dir / "representation_manifest.json"
    _write_json(representation_path, representation_manifest)

    mapped_latent = {}
    for block_name in ("gse128639", "gse254985", "hgsoc"):
        block = prepared[block_name]
        batches = block.get("batches", block.get("samples"))
        mapped_latent[block_name] = map_scvi_query_latent(
            block["counts"],
            batches=batches,
            gene_order=gene_order,
            reference_model_dir=reference_scvi_dir,
            output_model_dir=application_dir / f"scvi_{block_name}",
            max_epochs=0,
            seed=seed,
            progress=False,
        )

    predictions = {}
    ridge_predictions = {}
    cognate_predictions = {}
    for block_name in ("gse128639", "gse254985", "hgsoc"):
        block = prepared[block_name]
        predictions[block_name] = predict_protrans(
            model,
            block["expression"],
            mapped_latent[block_name],
            gene_embedding,
            protein_embedding,
            batch_size=int(args.batch_size),
            device=device,
        )
        ridge_predictions[block_name] = predict_ridge_baseline(
            ridge, block["expression"]
        )
        cognate_predictions[block_name] = cognate_rna_score(
            block["expression"],
            gene_order=gene_order,
            target="HLA_DR_COMPLEX",
        )

    hla_model_index = targets.index("HLA_DR_COMPLEX")
    independent_units = [
        (
            "GSE128639_all",
            "gse128639",
            np.ones(len(prepared["gse128639"]["cells"]), dtype=bool),
        )
    ]
    for condition in np.unique(prepared["gse254985"]["conditions"]):
        independent_units.append(
            (
                f"GSE254985_{condition}",
                "gse254985",
                np.asarray(prepared["gse254985"]["conditions"]) == condition,
            )
        )
    independent_rows = []
    independent_model_metrics = []
    independent_baseline_metrics = []
    for unit, block_name, mask in independent_units:
        block = prepared[block_name]
        block_targets = [str(value) for value in block["targets"]]
        truth_index = block_targets.index("HLA_DR_COMPLEX")
        truth = np.asarray(block["protein"])[mask, truth_index]
        method_predictions = {
            "scprotrans": predictions[block_name][mask, hla_model_index],
            "ridge": ridge_predictions[block_name][mask, hla_model_index],
            "cognate_rna": cognate_predictions[block_name][mask],
        }
        unit_metrics = {}
        for method, values in method_predictions.items():
            metrics = compute_prediction_metrics(truth, values)
            unit_metrics[method] = metrics
            independent_rows.append(
                {
                    "evaluation": "independent_transfer",
                    "unit": unit,
                    "target": "HLA_DR_COMPLEX",
                    "method": method,
                    "cell_count": int(np.sum(mask)),
                    **metrics,
                }
            )
        independent_model_metrics.append(unit_metrics["scprotrans"]["spearman"])
        independent_baseline_metrics.append(
            max(
                unit_metrics["ridge"]["spearman"],
                unit_metrics["cognate_rna"]["spearman"],
            )
        )
    _write_tsv(application_dir / "independent_transfer_metrics.tsv", independent_rows)

    ood_model = fit_knn_ood(
        reference_latent,
        reference_groups=reference["donors"],
        k=10,
        quantile=0.95,
    )
    hgsoc_scores, hgsoc_in_domain = score_knn_ood(
        ood_model, mapped_latent["hgsoc"]
    )
    hgsoc = prepared["hgsoc"]
    strict_mask = np.asarray(hgsoc["high_confidence"], dtype=bool)
    if not np.any(strict_mask):
        raise ValueError("The HGSOC prepared block contains no strict-EOC cells.")
    ood_summary = patient_equal_in_domain_fraction(
        hgsoc_in_domain[strict_mask], np.asarray(hgsoc["patients"])[strict_mask]
    )
    ood_path = application_dir / "knn_ood_model.npz"
    np.savez_compressed(
        ood_path,
        reference_latent=np.asarray(ood_model["reference_latent"]),
        k=np.asarray([ood_model["k"]]),
        quantile=np.asarray([ood_model["quantile"]]),
        threshold=np.asarray([ood_model["threshold"]]),
    )

    eligible_path = (
        Path(args.output_dir).parent
        / "tables"
        / "gse266577_patient_analysis_sets.tsv"
    )
    eligible_patients = read_eligible_patients(eligible_path)
    paired_mask = np.isin(hgsoc["patients"], sorted(eligible_patients))
    paired = compute_patient_paired_effect(
        predictions["hgsoc"][paired_mask, hla_model_index],
        np.asarray(hgsoc["patients"])[paired_mask],
        np.asarray(hgsoc["stages"])[paired_mask],
        pre_label="chemo-naive",
        post_label="IDS",
    )
    patient_rows = [
        {
            "patient_id": patient,
            "predicted_hla_dr_post_minus_pre": float(effect),
            "seed": seed,
        }
        for patient, effect in zip(paired["patients"], paired["patient_effects"])
    ]
    _write_tsv(application_dir / "hgsoc_patient_effects.tsv", patient_rows)
    hgsoc_rows = [
        {
            "cell": str(cell),
            "patient_id": str(patient),
            "sample_id": str(sample),
            "treatment_stage": str(stage),
            "strict_eoc": bool(strict),
            "predicted_hla_dr": float(prediction),
            "ood_score": float(score),
            "in_domain": bool(in_domain),
            "seed": seed,
        }
        for cell, patient, sample, stage, strict, prediction, score, in_domain in zip(
            hgsoc["cells"],
            hgsoc["patients"],
            hgsoc["samples"],
            hgsoc["stages"],
            hgsoc["high_confidence"],
            predictions["hgsoc"][:, hla_model_index],
            hgsoc_scores,
            hgsoc_in_domain,
        )
    ]
    _write_tsv(application_dir / "hgsoc_cell_predictions.tsv.gz", hgsoc_rows)

    bundle_path = application_dir / "scprotrans_model_bundle.pt"
    encoder_path = reference_scvi_dir / "model.pt"
    metadata = {
        "gene_order": gene_order,
        "protein_order": targets,
        "normalization": {
            "rna": "library_size_1e4_log1p",
            "adt": "targetwise_log1p_raw_counts",
        },
        "seed": seed,
        "model_config": model_config,
        "training_config": {
            "epochs": int(args.epochs),
            "completed_epochs": len(history),
            "batch_size": int(args.batch_size),
            "patience": int(args.patience),
        },
        "split_manifest": representation_manifest,
        "input_hashes": prepared["manifest"]["source_sha256"],
        "embedding_hashes": prepared["manifest"]["embedding_sha256"],
        "representation": {
            "encoder": os.fspath(encoder_path),
            "encoder_sha256": sha256_file(encoder_path),
            "query_adaptation_epochs": 0,
        },
        "ood_artifact": {
            "path": os.fspath(ood_path),
            "sha256": sha256_file(ood_path),
        },
    }
    save_model_bundle(bundle_path, model, metadata)
    restored_model, _ = load_model_bundle(bundle_path)
    roundtrip_expected = predict_protrans(
        model,
        reference["expression"][:32],
        reference_latent[:32],
        gene_embedding,
        protein_embedding,
        batch_size=16,
        device="cpu",
    )
    roundtrip_observed = predict_protrans(
        restored_model,
        reference["expression"][:32],
        reference_latent[:32],
        gene_embedding,
        protein_embedding,
        batch_size=16,
        device="cpu",
    )
    artifact_hashes = {
        os.fspath(bundle_path): sha256_file(bundle_path),
        os.fspath(encoder_path): sha256_file(encoder_path),
        os.fspath(ood_path): sha256_file(ood_path),
        os.fspath(representation_path): sha256_file(representation_path),
    }
    bundle_audit = {
        "roundtrip_match": bool(
            np.allclose(roundtrip_expected, roundtrip_observed, atol=1e-7)
        ),
        "artifact_hashes_match": bool(
            all(sha256_file(path) == digest for path, digest in artifact_hashes.items())
        ),
        "encoder_present": encoder_path.exists(),
        "split_manifest_present": representation_path.exists(),
        "artifact_sha256": artifact_hashes,
    }
    _write_json(application_dir / "bundle_audit.json", bundle_audit)
    summary = {
        "seed": seed,
        "independent_model_metrics": independent_model_metrics,
        "independent_baseline_metrics": independent_baseline_metrics,
        "independent_noninferiority": paired_noninferiority_summary(
            independent_model_metrics,
            independent_baseline_metrics,
            margin=0.05,
            iterations=1000,
            seed=seed,
        ),
        "hgsoc_ood_summary": ood_summary,
        "hgsoc_paired_effect": {
            "patient_count": paired["patient_count"],
            "mean_effect": paired["mean_effect"],
            "median_effect": paired["median_effect"],
            "patient_effect_sd": paired["patient_effect_sd"],
        },
        "bundle_audit": bundle_audit,
    }
    _write_json(application_dir / "application_summary.json", summary)
    return {
        **summary,
        "hgsoc_in_domain": hgsoc_in_domain[strict_mask],
        "hgsoc_patients": np.asarray(hgsoc["patients"])[strict_mask],
        "patient_effects": paired["patient_effects"],
    }


def run_pilot_benchmark(args: argparse.Namespace) -> dict[str, object]:
    prepared = prepare_multidomain_data(args)
    seed = int(args.seeds[0])
    config_tag = (
        f"h{args.model_hidden_dim}_sv{args.scvi_epochs}_pt{args.epochs}_"
        f"seed{seed}"
    )
    run_dir = Path(args.output_dir) / "pilot" / config_tag
    run_dir.mkdir(parents=True, exist_ok=True)
    donor = run_donor_holdout_benchmark(
        prepared,
        args,
        seed=seed,
        run_dir=run_dir,
    )
    donor_summary = donor["hla_dr_summary"]
    stop_loss_passed = bool(
        donor_summary["estimate"] >= 0.20
        and donor_summary["finite_donor_count"] >= 5
    )
    summary = {
        "run": "pilot",
        "seed": seed,
        "prepared_cache": os.fspath(prepared["cache_dir"]),
        "hla_dr_donor_holdout": donor_summary,
        "stop_loss_threshold": {
            "minimum_donor_equal_spearman": 0.20,
            "minimum_finite_donors": 5,
        },
        "stop_loss_passed": stop_loss_passed,
        "gate3_evaluated": False,
        "gate3_reason": (
            "Pilot passed; independent transfer and five-seed HGSOC stability remain required."
            if stop_loss_passed
            else "Pilot failed the preregistered stop-loss; formal scProTrans claims are not authorized."
        ),
    }
    _write_json(run_dir / "pilot_summary.json", summary)
    return summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run",
        choices=("contract", "pilot", "full"),
        default="contract",
        help="Run only contracts, one-seed stop-loss benchmark, or the frozen full benchmark.",
    )
    parser.add_argument(
        "--contract-only",
        action="store_true",
        help="Deprecated alias for --run contract.",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(
            "data/raw/scprotrans_reference"
        ),
    )
    parser.add_argument(
        "--hgsoc-root",
        type=Path,
        default=Path("data/raw/gse266577"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/scprotrans_hgsoc_v4/scprotrans_benchmark"),
    )
    parser.add_argument("--max-genes", type=int, default=512)
    parser.add_argument("--max-ref-cells-per-donor", type=int, default=1000)
    parser.add_argument("--max-hgsoc-cells", type=int, default=0)
    parser.add_argument("--scvi-latent", type=int, default=50)
    parser.add_argument("--scvi-epochs", type=int, default=15)
    parser.add_argument("--epochs", type=int, default=15)
    parser.add_argument("--patience", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--model-hidden-dim", type=int, default=100)
    parser.add_argument("--model-heads", type=int, default=5)
    parser.add_argument("--seeds", nargs="+", type=int, default=[260716])
    parser.add_argument("--device", choices=("auto", "cpu", "mps", "cuda"), default="auto")
    return parser.parse_args(argv)


def main() -> int:
    args = parse_args()
    if args.contract_only or args.run == "contract":
        model = ProTrans()
        print(f"ProTrans contract ready with {sum(p.numel() for p in model.parameters())} parameters")
        return 0
    if args.run == "pilot":
        print(json.dumps(run_pilot_benchmark(args), indent=2, sort_keys=True))
        return 0
    raise SystemExit(
        "The full five-seed runner is gated on a successful pilot benchmark."
    )


if __name__ == "__main__":
    raise SystemExit(main())
