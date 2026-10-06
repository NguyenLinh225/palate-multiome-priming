"""Cluster annotation from absolute marker expression.

Cross-cell-type z-scoring is deliberately avoided. A cluster whose markers are
all weakly expressed can still rank highest for some lineage after z-scoring,
which mislabels it; thresholds on mean log-normalized expression do not have
that failure mode.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd


def panel_scores(adata: ad.AnnData, panels: dict[str, list[str]],
                 groupby: str = "leiden") -> pd.DataFrame:
    """Mean log-normalized expression of each marker panel, per cluster."""
    source = adata.raw.to_adata() if adata.raw is not None else adata
    cols = {}
    for name, genes in panels.items():
        present = [g for g in genes if g in source.var_names]
        if not present:
            continue
        cols[name] = np.asarray(source[:, present].X.mean(1)).ravel()
    scores = pd.DataFrame(cols, index=adata.obs_names)
    scores[groupby] = adata.obs[groupby].to_numpy()
    return scores.groupby(groupby, observed=True).mean()


def detection_fraction(adata: ad.AnnData, gene: str, groupby: str = "leiden") -> pd.Series:
    """Fraction of cells per cluster with non-zero counts for one gene."""
    source = adata.raw.to_adata() if adata.raw is not None else adata
    if gene not in source.var_names:
        return pd.Series(0.0, index=adata.obs[groupby].cat.categories)
    v = np.asarray(source[:, [gene]].X.todense()).ravel() > 0
    return pd.Series(v, index=adata.obs_names).groupby(adata.obs[groupby].to_numpy()).mean()


def call_lineage(scores: pd.DataFrame, thresholds: dict) -> dict[str, str]:
    """Assign each cluster a major lineage; most specific marker wins."""
    order = ["Erythroid", "Endothelial", "Immune", "Neural/glia", "Myocyte"]
    calls = {}
    for cluster, row in scores.iterrows():
        label = "Other"
        for lineage in order:
            if lineage in row and row[lineage] >= thresholds[lineage]:
                label = lineage
                break
        else:
            if row.get("Epithelial", 0) >= thresholds["Epithelial"] and \
                    row.get("Mesenchymal", 0) < thresholds["Mesenchymal"]:
                label = "Epithelial"
            elif row.get("Mesenchymal", 0) >= thresholds["Mesenchymal"]:
                label = "Mesenchymal"
        calls[cluster] = label
    return calls


def call_mesenchymal_subtype(scores: pd.DataFrame, thresholds: dict,
                             acan_fraction: pd.Series | None = None) -> dict[str, str]:
    """Sub-label mesenchymal clusters as anterior, posterior, osteogenic or early.

    A chondrogenic call requires *Acan* detection above
    ``chondrogenic_min_acan_frac``; Sox9/Col2a1 alone are expressed in
    uncommitted palatal mesenchyme and would over-call chondrocytes.
    """
    calls = {}
    for cluster, row in scores.iterrows():
        ant, post = row.get("Anterior", 0.0), row.get("Posterior", 0.0)
        if row.get("Osteogenic", 0.0) >= thresholds["osteogenic_min"]:
            calls[cluster] = "Osteogenic"
        elif ant >= thresholds["axis_min"] and ant > thresholds["axis_fold_over_opposite"] * post:
            calls[cluster] = "Anterior"
        elif post >= thresholds["axis_min"] and post > thresholds["axis_fold_over_opposite"] * ant:
            calls[cluster] = "Posterior"
        elif (row.get("Chondrogenic", 0.0) >= thresholds["chondrogenic_min"]
              and acan_fraction is not None
              and acan_fraction.get(cluster, 0.0) >= thresholds["chondrogenic_min_acan_frac"]):
            calls[cluster] = "Chondrogenic"
        else:
            calls[cluster] = "Early/unpatterned mesenchyme"
    return calls


def annotate(adata: ad.AnnData, cfg: dict) -> ad.AnnData:
    """Write a ``celltype`` column combining lineage and mesenchymal subtype."""
    lineage = panel_scores(adata, cfg["markers"]["lineage"])
    major = call_lineage(lineage, cfg["lineage_thresholds"])

    subtype_scores = panel_scores(adata, cfg["markers"]["mesenchymal_subtype"])
    acan = detection_fraction(adata, "Acan")
    subtype = call_mesenchymal_subtype(subtype_scores, cfg["subtype_thresholds"], acan)

    labels = {c: (subtype[c] if major[c] == "Mesenchymal" else major[c]) for c in major}
    order = ["Anterior", "Posterior", "Osteogenic", "Chondrogenic",
             "Early/unpatterned mesenchyme", "Epithelial", "Neural/glia", "Myocyte",
             "Endothelial", "Immune", "Erythroid", "Other"]
    values = adata.obs["leiden"].map(labels)
    present = [o for o in order if o in set(values)]
    adata.obs["celltype"] = pd.Categorical(values, categories=present, ordered=True)
    return adata
