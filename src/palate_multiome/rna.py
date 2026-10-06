"""RNA normalization, HVG selection, embedding and clustering."""
from __future__ import annotations

import warnings

import anndata as ad
import numpy as np
import scanpy as sc

S_GENES = ["Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6",
           "Cdca7", "Dtl", "Prim1", "Uhrf1", "Cenpu", "Hells", "Rfc2", "Rpa2", "Nasp",
           "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2",
           "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", "Casp8ap2",
           "Usp1", "Clspn", "Pola1", "Chaf1b", "Brip1", "E2f8"]
G2M_GENES = ["Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80", "Cks2",
             "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf", "Tacc3", "Smc4", "Ccnb2", "Ckap2l",
             "Ckap2", "Aurkb", "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b",
             "Hjurp", "Cdca3", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", "Ncapd2",
             "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr", "Aurka", "Psrc1",
             "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa"]


def normalize(adata: ad.AnnData, cfg: dict) -> ad.AnnData:
    """Library-size normalize, log1p, and stash the result in ``.raw``."""
    sc.pp.normalize_total(adata, target_sum=cfg["rna"]["target_sum"])
    sc.pp.log1p(adata)
    adata.raw = adata
    return adata


MIN_CYCLE_GENES = 5


def score_cell_cycle(adata: ad.AnnData) -> ad.AnnData:
    """Attach S/G2M scores so cycling can be inspected as a confounder.

    Scores are diagnostic, not used downstream, so a gene panel that lacks the
    cell-cycle genes (a different reference, a targeted panel) yields NaN scores
    and phase ``"unscored"`` with a warning instead of aborting the pipeline.
    """
    s = [g for g in S_GENES if g in adata.var_names]
    g2m = [g for g in G2M_GENES if g in adata.var_names]
    if len(s) < MIN_CYCLE_GENES or len(g2m) < MIN_CYCLE_GENES:
        warnings.warn(f"only {len(s)} S and {len(g2m)} G2M genes present; "
                      "skipping cell-cycle scoring", stacklevel=2)
        adata.obs["S_score"] = np.nan
        adata.obs["G2M_score"] = np.nan
        adata.obs["phase"] = "unscored"
        return adata
    sc.tl.score_genes_cell_cycle(adata, s_genes=s, g2m_genes=g2m)
    return adata


def embed(adata: ad.AnnData, cfg: dict) -> ad.AnnData:
    """HVG -> scale -> PCA -> neighbors -> UMAP -> Leiden.

    HVGs are selected per library (``hvg_batch_key``) so genes driven by a
    single library are down-weighted without an explicit integration step.
    """
    rna, seed = cfg["rna"], cfg["seed"]
    sc.pp.highly_variable_genes(adata, n_top_genes=rna["n_top_genes"],
                                batch_key=rna["hvg_batch_key"])
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=rna["scale_max_value"])
    sc.tl.pca(adata, n_comps=rna["n_pcs"], svd_solver="arpack", random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=rna["neighbors_k"],
                    n_pcs=rna["neighbors_n_pcs"], random_state=seed)
    sc.tl.umap(adata, min_dist=rna["umap_min_dist"], random_state=seed)
    sc.tl.leiden(adata, resolution=rna["leiden_resolution"], key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False, random_state=seed)
    return adata
