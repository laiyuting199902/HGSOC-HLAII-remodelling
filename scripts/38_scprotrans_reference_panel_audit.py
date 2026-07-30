#!/usr/bin/env python3
"""Audit scProTrans reference panels and target-protein feasibility."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import struct
import subprocess
import zipfile
from collections import OrderedDict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np


TARGET_SPECS = OrderedDict(
    [
        (
            "HLA_DR_COMPLEX",
            {
                "display_name": "HLA-DR/MHC-II complex",
                "priority": 1,
                "adt_resolution": "heterodimer_or_complex_adt",
                "sequence_strategy": "mean_chain_embedding",
                "query_genes": ("HLA-DRA", "HLA-DRB1"),
                "query_uniprot_ids": ("P01903", "P01911"),
                "claim_boundary": (
                    "The ADT measures an HLA-DR complex epitope and cannot be "
                    "reported as HLA-DRA or HLA-DRB1 single-chain abundance."
                ),
            },
        ),
        (
            "HLA_II_PAN",
            {
                "display_name": "pan HLA-DR/DP/DQ",
                "priority": 2,
                "adt_resolution": "pan_mhcii_complex_adt",
                "sequence_strategy": "mean_chain_embedding",
                "query_genes": (
                    "HLA-DRA",
                    "HLA-DRB1",
                    "HLA-DPA1",
                    "HLA-DPB1",
                    "HLA-DQA1",
                    "HLA-DQB1",
                ),
                "query_uniprot_ids": (
                    "P01903",
                    "P01911",
                    "P20036",
                    "P04440",
                    "P01909",
                    "P01920",
                ),
                "claim_boundary": (
                    "The pan-HLA-II ADT does not resolve DR, DP, DQ, or either "
                    "chain separately."
                ),
            },
        ),
        (
            "HLA_I_PAN",
            {
                "display_name": "pan HLA class I",
                "priority": 3,
                "adt_resolution": "pan_hlai_complex_adt",
                "sequence_strategy": "mean_chain_embedding",
                "query_genes": ("HLA-A", "HLA-B", "HLA-C", "B2M"),
                "query_uniprot_ids": ("P04439", "P01889", "P10321", "P61769"),
                "claim_boundary": (
                    "The pan-HLA-I ADT measures a class-I complex epitope and "
                    "does not resolve HLA-A, HLA-B, HLA-C, or B2M."
                ),
            },
        ),
        (
            "B2M",
            {
                "display_name": "B2M",
                "priority": 4,
                "adt_resolution": "single_protein_adt",
                "sequence_strategy": "single_protein_embedding",
                "query_genes": ("B2M",),
                "query_uniprot_ids": ("P61769",),
                "claim_boundary": "Only a B2M-specific ADT supports a B2M protein claim.",
            },
        ),
        (
            "PD_L1",
            {
                "display_name": "PD-L1/CD274",
                "priority": 5,
                "adt_resolution": "single_protein_adt",
                "sequence_strategy": "single_protein_embedding",
                "query_genes": ("CD274",),
                "query_uniprot_ids": ("Q9NZQ7",),
                "claim_boundary": "Prediction remains ADT-like abundance, not measured protein.",
            },
        ),
        (
            "EPCAM",
            {
                "display_name": "EpCAM/CD326",
                "priority": 6,
                "adt_resolution": "single_protein_adt",
                "sequence_strategy": "single_protein_embedding",
                "query_genes": ("EPCAM",),
                "query_uniprot_ids": ("P16422",),
                "claim_boundary": "Prediction remains ADT-like abundance, not measured protein.",
            },
        ),
        (
            "CD47",
            {
                "display_name": "CD47",
                "priority": 7,
                "adt_resolution": "single_protein_adt",
                "sequence_strategy": "single_protein_embedding",
                "query_genes": ("CD47",),
                "query_uniprot_ids": ("Q08722",),
                "claim_boundary": "Prediction remains ADT-like abundance, not measured protein.",
            },
        ),
        (
            "CD74",
            {
                "display_name": "CD74/invariant chain",
                "priority": 8,
                "adt_resolution": "single_protein_adt",
                "sequence_strategy": "single_protein_embedding",
                "query_genes": ("CD74",),
                "query_uniprot_ids": ("P04233",),
                "claim_boundary": (
                    "CD74 enters protein inference only after a real CD74 ADT "
                    "benchmark; RNA expression alone is insufficient."
                ),
            },
        ),
    ]
)


def _compact_label(label: str) -> str:
    return "".join(ch for ch in str(label).upper() if ch.isalnum())


def canonical_target(label: str) -> str | None:
    """Map heterogeneous antibody labels to evidence-aware target classes."""
    raw = str(label).strip()
    compact = _compact_label(raw)
    if not compact or "CTRL" in compact or "CONTROL" in compact or "IGG" in compact:
        return None
    if "HLADRDPDQ" in compact:
        return "HLA_II_PAN"
    if "HLAABC" in compact or "HLAABCPAN" in compact:
        return "HLA_I_PAN"
    if "MHCII" in compact or compact.startswith("HLADR"):
        return "HLA_DR_COMPLEX"
    if compact.startswith("CD74"):
        return "CD74"
    if compact.startswith("B2M") or "BETA2MICROGLOBULIN" in compact:
        return "B2M"
    if "PDL1" in compact or compact.startswith("CD274"):
        return "PD_L1"
    if "EPCAM" in compact or compact.startswith("CD326"):
        return "EPCAM"
    if compact.startswith("CD47"):
        return "CD47"
    return None


def _open_text(path: os.PathLike[str] | str):
    path = os.fspath(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="", encoding="utf-8")


def read_row_feature_labels(
    path: os.PathLike[str] | str,
    delimiter: str,
) -> tuple[list[str], int]:
    """Stream feature-by-cell text matrices without loading numeric values."""
    labels: list[str] = []
    with _open_text(path) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader)
        cell_count = len(header) - 1 if header and header[0] == "" else len(header)
        for row in reader:
            if row and row[0].strip():
                labels.append(row[0].strip())
    return labels, cell_count


def read_row_feature_matrix_index(
    path: os.PathLike[str] | str,
    delimiter: str,
) -> tuple[list[str], list[str]]:
    """Read feature and cell identifiers from a feature-by-cell text matrix."""
    labels: list[str] = []
    with _open_text(path) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader)
        cells = header[1:] if header and header[0] == "" else header
        for row in reader:
            if row and row[0].strip():
                labels.append(row[0].strip())
    return labels, [cell.strip() for cell in cells if cell.strip()]


def exclude_metadata_sidecars(paths: Iterable[os.PathLike[str] | str]) -> list[str]:
    """Remove Finder metadata files that appear beside downloaded GEO files."""
    return [
        os.fspath(path)
        for path in paths
        if not Path(path).name.startswith("._") and Path(path).name != ".DS_Store"
    ]


def rscript_available() -> bool:
    return shutil.which("Rscript") is not None


def read_double_gzip_rds_index(
    path: os.PathLike[str] | str,
) -> tuple[int, int, set[str]]:
    """Read dimensions and cell names from a doubly gzip-wrapped matrix RDS."""
    if not rscript_available():
        raise RuntimeError("Rscript is required to audit doubly compressed RDS matrices.")
    code = """
args <- commandArgs(TRUE)
con <- gzcon(gzfile(args[[1]], "rb"))
x <- readRDS(con)
close(con)
d <- dim(x)
if (length(d) != 2) stop("RDS object must be a two-dimensional matrix")
cat(d[[1]], d[[2]], sep="\\t")
cat("\\n")
cn <- colnames(x)
if (!is.null(cn)) writeLines(cn)
"""
    result = subprocess.run(
        ["Rscript", "-e", code, os.fspath(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    lines = result.stdout.splitlines()
    if not lines or "\t" not in lines[0]:
        raise ValueError(f"RDS audit did not return dimensions: {path}")
    rows_text, columns_text = lines[0].split("\t", 1)
    return int(rows_text), int(columns_text), set(lines[1:])


def summarize_cell_alignment(
    rna_cells: set[str],
    adt_cells: set[str],
    csp_cells: set[str],
    hto_cells: set[str],
) -> dict[str, int]:
    """Summarize complete-case cell overlap across CITE-seq modalities."""
    complete = rna_cells.intersection(adt_cells, csp_cells, hto_cells)
    return {
        "complete_cell_count": len(complete),
        "rna_only_count": len(rna_cells - complete),
        "adt_only_count": len(adt_cells - complete),
        "csp_only_count": len(csp_cells - complete),
        "hto_only_count": len(hto_cells - complete),
    }


def read_10x_feature_labels(path: os.PathLike[str] | str) -> list[str]:
    labels: list[str] = []
    with _open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if row:
                labels.append((row[1] if len(row) > 1 else row[0]).strip())
    return labels


def summarize_panel(labels: Sequence[str]) -> dict[str, object]:
    mapped = [canonical_target(label) for label in labels]
    measured = {target for target in mapped if target is not None}
    return {
        "adt_count": len(labels),
        "target_adt_count": sum(target is not None for target in mapped),
        "measured_targets": measured,
        "has_hla_dr_mhcii": bool(
            measured.intersection({"HLA_DR_COMPLEX", "HLA_II_PAN"})
        ),
        "has_hla_i_b2m": bool(measured.intersection({"HLA_I_PAN", "B2M"})),
        "has_pd_l1": "PD_L1" in measured,
        "has_epcam": "EPCAM" in measured,
        "has_cd47": "CD47" in measured,
        "has_cd74": "CD74" in measured,
    }


def _decode_strings(values: np.ndarray) -> set[str]:
    result = set()
    for value in np.asarray(values).reshape(-1):
        if isinstance(value, bytes):
            result.add(value.decode("utf-8"))
        else:
            result.add(str(value))
    return result


def _stored_npy_member_layout(
    path: os.PathLike[str] | str,
    member_name: str,
) -> tuple[tuple[int, ...], bool, np.dtype, int]:
    """Locate an uncompressed NPY member inside NPZ without materializing it."""
    from numpy.lib import format as npy_format

    path = os.fspath(path)
    with zipfile.ZipFile(path) as archive:
        info = archive.getinfo(member_name)
        if info.compress_type != zipfile.ZIP_STORED:
            raise ValueError(f"NPZ member must be stored for mmap: {member_name}")
        header_offset = info.header_offset

    with open(path, "rb") as handle:
        handle.seek(header_offset)
        local_header = handle.read(30)
        if len(local_header) != 30 or local_header[:4] != b"PK\x03\x04":
            raise ValueError(f"Invalid local ZIP header for {member_name}")
        filename_length, extra_length = struct.unpack("<HH", local_header[26:30])
        npy_offset = header_offset + 30 + filename_length + extra_length
        handle.seek(npy_offset)
        version = npy_format.read_magic(handle)
        if version == (1, 0):
            shape, fortran_order, dtype = npy_format.read_array_header_1_0(handle)
        elif version == (2, 0):
            shape, fortran_order, dtype = npy_format.read_array_header_2_0(handle)
        else:
            shape, fortran_order, dtype = npy_format.read_array_header_2_0(handle)
        data_offset = handle.tell()
    return tuple(shape), bool(fortran_order), np.dtype(dtype), data_offset


def _read_fixed_unicode_names(
    path: os.PathLike[str] | str,
    shape: tuple[int, ...],
    dtype: np.dtype,
    data_offset: int,
) -> set[str]:
    if len(shape) != 1 or dtype.kind != "U":
        raise ValueError("Gene member must be a one-dimensional fixed Unicode array.")
    raw = np.memmap(
        path,
        mode="r",
        dtype=np.uint8,
        offset=data_offset,
        shape=(shape[0], dtype.itemsize),
        order="C",
    )
    byte_order = "utf-32-be" if dtype.byteorder == ">" else "utf-32-le"
    probe_bytes = min(dtype.itemsize, 1024)
    names: set[str] = set()
    for index in range(shape[0]):
        probe = bytes(raw[index, :probe_bytes])
        text = probe.decode(byte_order, errors="ignore").split("\x00", 1)[0]
        if text:
            names.add(text)
    return names


def read_gene_embedding_index(path: os.PathLike[str] | str) -> tuple[set[str], int]:
    gene_shape, gene_fortran, gene_dtype, gene_offset = _stored_npy_member_layout(
        path, "gene.npy"
    )
    embedding_shape, _, _, _ = _stored_npy_member_layout(path, "embedding.npy")
    if gene_fortran or len(embedding_shape) != 2 or embedding_shape[0] != gene_shape[0]:
        raise ValueError("Gene embedding matrix is not gene-by-dimension.")
    genes = _read_fixed_unicode_names(path, gene_shape, gene_dtype, gene_offset)
    return genes, int(embedding_shape[1])


def read_protein_embedding_index(path: os.PathLike[str] | str) -> tuple[set[str], int]:
    import h5py

    with h5py.File(path, "r") as handle:
        entries = set(handle.keys())
        if not entries:
            raise ValueError("Protein embedding file is empty.")
        dimensions = {int(np.asarray(handle[key]).reshape(-1).shape[0]) for key in entries}
    if len(dimensions) != 1:
        raise ValueError("Protein embeddings have inconsistent dimensions.")
    return entries, dimensions.pop()


def build_target_coverage_rows(
    measured_targets: Mapping[str, set[str]],
    protein_entries: set[str],
    gene_names: set[str],
    hgsoc_gene_names: set[str] | None = None,
    dataset_donor_counts: Mapping[str, int] | None = None,
) -> list[dict[str, object]]:
    rows = []
    for target_id, spec in TARGET_SPECS.items():
        measured_arms = sorted(
            arm for arm, targets in measured_targets.items() if target_id in targets
        )
        protein_complete = all(
            entry in protein_entries for entry in spec["query_uniprot_ids"]
        )
        gene_complete = all(gene in gene_names for gene in spec["query_genes"])
        hgsoc_complete = (
            all(gene in hgsoc_gene_names for gene in spec["query_genes"])
            if hgsoc_gene_names is not None
            else False
        )
        measured = bool(measured_arms)
        benchmark_eligible = measured and protein_complete and gene_complete
        donor_eligible = benchmark_eligible and (
            any(dataset_donor_counts.get(arm, 0) >= 3 for arm in measured_arms)
            if dataset_donor_counts is not None
            else True
        )
        rows.append(
            {
                "target_id": target_id,
                "display_name": spec["display_name"],
                "priority": spec["priority"],
                "evidence_class": (
                    "measured-protein" if measured else "unvalidated-extrapolation"
                ),
                "measured_dataset_count": len(measured_arms),
                "measured_dataset_arms": ";".join(measured_arms),
                "adt_resolution": spec["adt_resolution"],
                "sequence_strategy": spec["sequence_strategy"],
                "query_genes": ";".join(spec["query_genes"]),
                "query_uniprot_ids": ";".join(spec["query_uniprot_ids"]),
                "protein_embedding_complete": protein_complete,
                "gene_embedding_complete": gene_complete,
                "hgsoc_transcript_complete": hgsoc_complete,
                "donor_holdout_eligible": donor_eligible,
                "protein_holdout_eligible": benchmark_eligible,
                "zero_shot_status": (
                    "withheld-protein-zero-shot-eligible"
                    if benchmark_eligible
                    else "not-eligible"
                ),
                "hgsoc_inference_authorized": False,
                "claim_boundary": spec["claim_boundary"],
            }
        )
    return rows


def task5_gate_passes(
    measured_targets: Mapping[str, set[str]],
    protein_entries: set[str],
    gene_names: set[str],
) -> bool:
    rows = build_target_coverage_rows(measured_targets, protein_entries, gene_names)
    return any(
        row["target_id"] in {"HLA_DR_COMPLEX", "HLA_II_PAN"}
        and row["protein_holdout_eligible"]
        for row in rows
    )


def read_gene_list(path: os.PathLike[str] | str) -> set[str]:
    genes: set[str] = set()
    with _open_text(path) as handle:
        for line in handle:
            gene = line.rstrip("\r\n").split("\t", 1)[0]
            if gene:
                genes.add(gene)
    return genes


def sha256_file(path: os.PathLike[str] | str, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def count_csv_metadata(
    path: os.PathLike[str] | str,
    donor_column: str,
    batch_column: str,
) -> tuple[int, int, int]:
    cells = 0
    donors: set[str] = set()
    batches: set[str] = set()
    with _open_text(path) as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            cells += 1
            if row.get(donor_column):
                donors.add(row[donor_column])
            if row.get(batch_column):
                batches.add(row[batch_column])
    return cells, len(donors), len(batches)


def read_h5ad_panel(path: os.PathLike[str] | str) -> tuple[list[str], int, int, int]:
    """Read H5AD axes directly with h5py; numeric matrices stay on disk."""
    import h5py

    def strings(dataset) -> list[str]:
        values = np.asarray(dataset).reshape(-1)
        return [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values]

    def column_strings(group, column: str) -> list[str]:
        node = group[column]
        if hasattr(node, "keys") and "categories" in node:
            categories = strings(node["categories"])
            codes = np.asarray(node["codes"]).reshape(-1)
            return [categories[int(code)] for code in codes if int(code) >= 0]
        categories_group = group.get("__categories")
        if categories_group is not None and column in categories_group:
            categories = strings(categories_group[column])
            codes = np.asarray(node).reshape(-1)
            return [categories[int(code)] for code in codes if int(code) >= 0]
        return strings(node)

    with h5py.File(path, "r") as handle:
        obs_names = strings(handle["obs"]["_index"])
        var_names = strings(handle["var"]["_index"])
        if "feature_types" in handle["var"]:
            feature_types = column_strings(handle["var"], "feature_types")
            labels = [
                name
                for name, feature_type in zip(var_names, feature_types)
                if feature_type.upper() in {"ADT", "ANTIBODY CAPTURE", "PROTEIN"}
            ]
        else:
            labels = [name for name in var_names if canonical_target(name) is not None]

        donor_count = 0
        batch_count = 0
        for column in ("DonorNumber", "DonorID", "donor", "donor_id"):
            if column in handle["obs"]:
                donor_count = len(set(column_strings(handle["obs"], column)))
                break
        for column in ("batch", "Batch", "Samplename"):
            if column in handle["obs"]:
                batch_count = len(set(column_strings(handle["obs"], column)))
                break
    return labels, len(obs_names), donor_count, batch_count


def audit_gse254985_arm(
    dataset_arm: str,
    accession: str,
    treatment: str,
    rna_path: os.PathLike[str] | str,
    adt_path: os.PathLike[str] | str,
    csp_path: os.PathLike[str] | str,
    hto_path: os.PathLike[str] | str,
) -> tuple[dict[str, object], set[str]]:
    """Audit one pancreatic-islet epithelial RNA/ADT/CSP/HTO arm."""
    rna_rows, rna_columns, rna_cells = read_double_gzip_rds_index(rna_path)
    adt_labels, adt_cells_list = read_row_feature_matrix_index(adt_path, delimiter="\t")
    csp_labels, csp_cells_list = read_row_feature_matrix_index(csp_path, delimiter="\t")
    _, hto_cells_list = read_row_feature_matrix_index(hto_path, delimiter="\t")
    adt_cells = set(adt_cells_list)
    csp_cells = set(csp_cells_list)
    hto_cells = set(hto_cells_list)
    alignment = summarize_cell_alignment(rna_cells, adt_cells, csp_cells, hto_cells)
    summary = summarize_panel(adt_labels + csp_labels)
    caveat = (
        "Human epithelial cytokine-induction bridge with measured HLA-DR ADT; "
        "it supports cross-domain transfer but is not an ovarian tumor cohort."
    )
    if treatment == "mixed_untreated_and_cytokine_treated":
        caveat += " Treatment assignment requires HTO demultiplexing before model evaluation."
    row = {
        "dataset_arm": dataset_arm,
        "accession": accession,
        "biological_context": "human_pancreatic_islet_epithelial",
        "treatment": treatment,
        "donor_count": 1,
        "batch_count": 1,
        "cell_count": alignment["complete_cell_count"],
        "rna_feature_count": rna_rows,
        "rna_cell_count": rna_columns,
        "complete_cell_count": alignment["complete_cell_count"],
        "rna_only_count": alignment["rna_only_count"],
        "adt_only_count": alignment["adt_only_count"],
        "csp_only_count": alignment["csp_only_count"],
        "hto_only_count": alignment["hto_only_count"],
        "adt_count": summary["adt_count"],
        "target_adt_count": summary["target_adt_count"],
        "has_hla_dr_mhcii": summary["has_hla_dr_mhcii"],
        "has_hla_i_b2m": summary["has_hla_i_b2m"],
        "has_pd_l1": summary["has_pd_l1"],
        "has_epcam": summary["has_epcam"],
        "has_cd47": summary["has_cd47"],
        "has_cd74": summary["has_cd74"],
        "measured_targets": ";".join(sorted(summary["measured_targets"])),
        "reference_role": "epithelial_cytokine_induction_bridge",
        "source_status": "downloaded_integrity_checked_and_audited",
        "source_file": ";".join(
            os.fspath(path) for path in (rna_path, adt_path, csp_path, hto_path)
        ),
        "caveat": caveat,
    }
    return row, set(summary["measured_targets"])


def audit_default_gse254985(
    data_root: Path,
) -> tuple[list[dict[str, object]], dict[str, set[str]]]:
    extracted = data_root / "GSE254985" / "extracted"
    if not extracted.exists():
        return [], {}
    available = {
        Path(path).name: Path(path)
        for path in exclude_metadata_sidecars(extracted.iterdir())
        if Path(path).is_file()
    }
    configs = [
        {
            "dataset_arm": "GSE254985_HP20276_untreated",
            "accession": "GSE254985/GSM8061741,GSM8061744,GSM8061747",
            "treatment": "untreated",
            "rna": "GSM8061741_20276-no_umiCleanMerged.rds.gz",
            "adt": "GSM8061744_20276-no_ADT_counts_n1000.txt.gz",
            "csp": "GSM8061744_20276-no_CSP_counts_n1000.txt.gz",
            "hto": "GSM8061747_20276-no_HTO_counts_n1000.txt.gz",
        },
        {
            "dataset_arm": "GSE254985_HP20276_cytokine_treated",
            "accession": "GSE254985/GSM8061742,GSM8061745,GSM8061748",
            "treatment": "48h_IFNG_IL1B_TNFA",
            "rna": "GSM8061742_20276plus_umiCleanMerged.rds.gz",
            "adt": "GSM8061745_20276plus_ADT_counts_n1000.txt.gz",
            "csp": "GSM8061745_20276plus_CSP_counts_n1000.txt.gz",
            "hto": "GSM8061748_20276plus_HTO_counts_n1000.txt.gz",
        },
        {
            "dataset_arm": "GSE254985_HP21024_mixed",
            "accession": "GSE254985/GSM8061743,GSM8061746,GSM8061749",
            "treatment": "mixed_untreated_and_cytokine_treated",
            "rna": "GSM8061743_21024_umiCleanMerged.rds.gz",
            "adt": "GSM8061746_21024_ADT_counts_n1000.txt.gz",
            "csp": "GSM8061746_21024_CSP_counts_n1000.txt.gz",
            "hto": "GSM8061749_21024_HTO_counts_n1000.txt.gz",
        },
    ]
    rows: list[dict[str, object]] = []
    measured_targets: dict[str, set[str]] = {}
    for config in configs:
        required = [config[key] for key in ("rna", "adt", "csp", "hto")]
        missing = [name for name in required if name not in available]
        if missing:
            raise FileNotFoundError(
                f"Incomplete GSE254985 arm {config['dataset_arm']}: {', '.join(missing)}"
            )
        row, targets = audit_gse254985_arm(
            dataset_arm=str(config["dataset_arm"]),
            accession=str(config["accession"]),
            treatment=str(config["treatment"]),
            rna_path=available[str(config["rna"])],
            adt_path=available[str(config["adt"])],
            csp_path=available[str(config["csp"])],
            hto_path=available[str(config["hto"])],
        )
        rows.append(row)
        measured_targets[str(config["dataset_arm"])] = targets
    return rows, measured_targets


def _bool_text(value: object) -> str:
    return "TRUE" if bool(value) else "FALSE"


def _write_tsv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(dict.fromkeys(key for row in rows for key in row.keys()))
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: _bool_text(value) if isinstance(value, bool) else value
                    for key, value in row.items()
                }
            )


def audit_default_datasets(data_root: Path) -> tuple[list[dict[str, object]], dict[str, set[str]]]:
    configs = [
        {
            "dataset_arm": "GSE164378_3P",
            "accession": "GSE164378/GSM5008738",
            "context": "PBMC",
            "role": "donor_holdout_training_reference",
            "reader": "10x",
            "panel": data_root / "GSE164378" / "GSM5008738_ADT_3P-features.tsv.gz",
            "metadata": data_root / "GSE164378" / "GSE164378_sc.meta.data_3P.csv.gz",
            "donors": 8,
            "batches": 2,
            "caveat": "Immune reference; HLA-DR is a complex ADT and no epithelial tumor domain is represented.",
        },
        {
            "dataset_arm": "GSE164378_5P",
            "accession": "GSE164378/GSM5008741",
            "context": "PBMC",
            "role": "donor_holdout_training_reference",
            "reader": "10x",
            "panel": data_root / "GSE164378" / "GSM5008741_ADT_5P-features.tsv.gz",
            "metadata": data_root / "GSE164378" / "GSE164378_sc.meta.data_5P.csv.gz",
            "donors": 8,
            "batches": 2,
            "caveat": "EpCAM is present in the antibody panel, but the assayed cells are PBMC rather than epithelial tumor cells.",
        },
        {
            "dataset_arm": "GSE128639_BMMC",
            "accession": "GSE128639/GSM3681519",
            "context": "bone_marrow_mononuclear_cells",
            "role": "independent_hla_dr_reference",
            "reader": "rows_tsv",
            "panel": data_root / "GSE128639" / "GSM3681519_MNC_ADT_counts.tsv.gz",
            "donors": 1,
            "batches": 1,
            "caveat": "Single source specimen; suitable for measured-target transfer, not donor holdout.",
        },
        {
            "dataset_arm": "GSE200417_CITE",
            "accession": "GSE200417/GSM6032898",
            "context": "activated_peripheral_blood_t_cells",
            "role": "cross_technology_reference",
            "reader": "rows_csv",
            "panel": data_root / "GSE200417" / "GSM6032898_CITE_ADT.csv.gz",
            "donors": 0,
            "batches": 1,
            "caveat": "Activated T-cell study; donor labels are not encoded in the downloaded ADT matrix.",
        },
        {
            "dataset_arm": "GSE200417_DOGMA",
            "accession": "GSE200417/GSM6032894",
            "context": "activated_peripheral_blood_t_cells",
            "role": "cross_technology_transfer",
            "reader": "rows_csv",
            "panel": data_root / "GSE200417" / "GSM6032894_DOGMA_ADT.csv.gz",
            "donors": 0,
            "batches": 1,
            "caveat": "Technology-transfer arm from the same study; not an independent donor cohort.",
        },
        {
            "dataset_arm": "GSE194122_CITE",
            "accession": "GSE194122",
            "context": "bone_marrow_mononuclear_cells",
            "role": "independent_multi_donor_transfer",
            "reader": "h5ad",
            "panel": data_root
            / "GSE194122"
            / "GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad",
            "donors": 6,
            "batches": 4,
            "caveat": "Multi-donor immune reference; target HGSOC epithelial cells still require an explicit OOD gate.",
        },
    ]

    rows: list[dict[str, object]] = []
    measured_targets: dict[str, set[str]] = {}
    for config in configs:
        path = Path(config["panel"])
        if not path.exists():
            labels: list[str] = []
            cell_count = 0
            source_status = "missing_or_pending"
            donor_count = int(config["donors"])
            batch_count = int(config["batches"])
        else:
            if config["reader"] == "10x":
                labels = read_10x_feature_labels(path)
                metadata = Path(config["metadata"])
                if metadata.exists():
                    cell_count, donor_count, batch_count = count_csv_metadata(
                        metadata, donor_column="donor", batch_column="Batch"
                    )
                else:
                    cell_count = 0
                    donor_count = int(config["donors"])
                    batch_count = int(config["batches"])
            elif config["reader"] == "rows_tsv":
                labels, cell_count = read_row_feature_labels(path, delimiter="\t")
                donor_count = int(config["donors"])
                batch_count = int(config["batches"])
            elif config["reader"] == "rows_csv":
                labels, cell_count = read_row_feature_labels(path, delimiter=",")
                donor_count = int(config["donors"])
                batch_count = int(config["batches"])
            elif config["reader"] == "h5ad":
                labels, cell_count, donor_count_h5, batch_count_h5 = read_h5ad_panel(path)
                donor_count = donor_count_h5 or int(config["donors"])
                batch_count = batch_count_h5 or int(config["batches"])
            else:
                raise ValueError(f"Unknown reader: {config['reader']}")
            source_status = "downloaded_and_audited"

        summary = summarize_panel(labels)
        measured_targets[config["dataset_arm"]] = set(summary["measured_targets"])
        rows.append(
            {
                "dataset_arm": config["dataset_arm"],
                "accession": config["accession"],
                "biological_context": config["context"],
                "donor_count": donor_count,
                "batch_count": batch_count,
                "cell_count": cell_count,
                "adt_count": summary["adt_count"],
                "target_adt_count": summary["target_adt_count"],
                "has_hla_dr_mhcii": summary["has_hla_dr_mhcii"],
                "has_hla_i_b2m": summary["has_hla_i_b2m"],
                "has_pd_l1": summary["has_pd_l1"],
                "has_epcam": summary["has_epcam"],
                "has_cd47": summary["has_cd47"],
                "has_cd74": summary["has_cd74"],
                "measured_targets": ";".join(sorted(summary["measured_targets"])),
                "reference_role": config["role"],
                "source_status": source_status,
                "source_file": str(path),
                "caveat": config["caveat"],
            }
        )
    epithelial_rows, epithelial_targets = audit_default_gse254985(data_root)
    rows.extend(epithelial_rows)
    measured_targets.update(epithelial_targets)
    return rows, measured_targets


def parse_args() -> argparse.Namespace:
    default_root = Path(
        "data/raw/scprotrans_reference"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=default_root)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/scprotrans_hgsoc_v4/tables"),
    )
    parser.add_argument("--summary-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    protein_path = args.data_root / "embeddings" / "per-protein.h5"
    gene_path = args.data_root / "embeddings" / "dna2vec_1w.npz"
    if not protein_path.exists() or not gene_path.exists():
        missing = [str(path) for path in (protein_path, gene_path) if not path.exists()]
        raise FileNotFoundError("Missing required official embedding(s): " + ", ".join(missing))

    dataset_rows, measured_targets = audit_default_datasets(args.data_root)
    protein_entries, protein_dimension = read_protein_embedding_index(protein_path)
    gene_names, gene_dimension = read_gene_embedding_index(gene_path)
    hgsoc_feature_path = (
        args.data_root.parent
        / "gse266577"
        / "GSE266577_seurat_features.txt.gz"
    )
    if not hgsoc_feature_path.exists():
        raise FileNotFoundError(f"Missing HGSOC feature list: {hgsoc_feature_path}")
    hgsoc_gene_names = read_gene_list(hgsoc_feature_path)
    dataset_donor_counts = {
        str(row["dataset_arm"]): int(row["donor_count"]) for row in dataset_rows
    }
    target_rows = build_target_coverage_rows(
        measured_targets,
        protein_entries,
        gene_names,
        hgsoc_gene_names=hgsoc_gene_names,
        dataset_donor_counts=dataset_donor_counts,
    )
    gate_passed = task5_gate_passes(measured_targets, protein_entries, gene_names)
    protein_sha256 = sha256_file(protein_path)
    gene_sha256 = sha256_file(gene_path)

    for row in target_rows:
        row["protein_embedding_dimension"] = protein_dimension
        row["gene_embedding_dimension"] = gene_dimension
        row["task5_hla_ii_gate_passed"] = gate_passed
        row["protein_embedding_sha256"] = protein_sha256
        row["gene_embedding_sha256"] = gene_sha256

    dataset_path = args.output_dir / "scprotrans_reference_dataset_audit.tsv"
    target_path = args.output_dir / "scprotrans_target_protein_coverage.tsv"
    _write_tsv(dataset_path, dataset_rows)
    _write_tsv(target_path, target_rows)

    summary = {
        "dataset_table": str(dataset_path),
        "target_table": str(target_path),
        "dataset_arms": len(dataset_rows),
        "protein_embedding_entries": len(protein_entries),
        "protein_embedding_dimension": protein_dimension,
        "gene_embedding_entries": len(gene_names),
        "gene_embedding_dimension": gene_dimension,
        "hgsoc_gene_entries": len(hgsoc_gene_names),
        "hgsoc_gene_embedding_overlap": len(hgsoc_gene_names.intersection(gene_names)),
        "hgsoc_gene_embedding_overlap_fraction": (
            len(hgsoc_gene_names.intersection(gene_names)) / len(hgsoc_gene_names)
        ),
        "task5_hla_ii_gate_passed": gate_passed,
        "authorized_hgsoc_targets": [],
    }
    if args.summary_json:
        print(json.dumps(summary, indent=2, ensure_ascii=False))
    else:
        print(
            "Task 5 HLA-II feasibility gate: "
            + ("PASS" if gate_passed else "FAIL")
            + "; HGSOC inference remains unauthorized until Gate 3."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
