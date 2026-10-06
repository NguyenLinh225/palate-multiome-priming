"""Per-cell QC metrics and the published filter gates."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from .io import _read_h5, cell_id

COUNT_METRICS = ["nCount_RNA", "nFeature_RNA", "percent_mt", "nCount_ATAC"]


def compute_metrics(path, sample: str, stage: str, mt_prefix: str = "mt-") -> pd.DataFrame:
    """RNA and ATAC count metrics for every barcode in one library."""
    X, feature_type, names, _, barcodes = _read_h5(path)
    X = X.tocsr()
    g = feature_type == b"Gene Expression"
    p = feature_type == b"Peaks"
    Xg, Xp = X[g], X[p]

    n_rna = np.asarray(Xg.sum(0)).ravel()
    mt = np.char.startswith(names[g], mt_prefix)
    return pd.DataFrame({
        "cell_id": [cell_id(b, sample) for b in barcodes],
        "sample": sample,
        "stage": stage,
        "nCount_RNA": n_rna,
        "nFeature_RNA": np.asarray((Xg > 0).sum(0)).ravel(),
        "percent_mt": 100 * np.asarray(Xg[mt].sum(0)).ravel() / np.maximum(n_rna, 1),
        "nCount_ATAC": np.asarray(Xp.sum(0)).ravel(),
    })


def count_gate(metrics: pd.DataFrame, qc: dict) -> pd.Series:
    """Boolean pass/fail for the count-based gates.

    TSS enrichment and nucleosome signal are *not* applied here: both need the
    fragment files. Expect ~5% more cells than the published post-QC count.
    """
    return ((metrics["nCount_RNA"] < qc["max_count_rna"])
            & (metrics["nFeature_RNA"] < qc["max_feature_rna"])
            & (metrics["percent_mt"] < qc["max_percent_mt"])
            & (metrics["nCount_ATAC"] > qc["min_count_atac"])
            & (metrics["nCount_ATAC"] < qc["max_count_atac"]))


def gate_summary(metrics: pd.DataFrame, qc: dict) -> pd.DataFrame:
    """Per-library pass counts, for the feasibility table."""
    m = metrics.assign(passes=count_gate(metrics, qc))
    out = (m.groupby(["stage", "sample"], observed=True)
             .agg(input_cells=("passes", "size"), pass_count_qc=("passes", "sum"))
             .reset_index())
    out["pct_pass"] = (100 * out["pass_count_qc"] / out["input_cells"]).round(1)
    return out


def apply_gate(adata: ad.AnnData, metrics: pd.DataFrame, qc: dict) -> ad.AnnData:
    """Subset an AnnData to passing cells and attach the metrics as obs columns."""
    idx = metrics.set_index("cell_id")
    keep = count_gate(idx, qc).reindex(adata.obs_names).fillna(False).to_numpy()
    out = adata[keep].copy()
    for c in COUNT_METRICS:
        out.obs[c] = idx.loc[out.obs_names, c].to_numpy()
    return out
