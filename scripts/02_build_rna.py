#!/usr/bin/env python
"""Step 2: build the annotated RNA object.

Loads all libraries, applies the count-based QC gates, normalizes, embeds and
annotates. Writes data/processed/palate_rna.h5ad and the cell-type tables.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from palate_multiome.annotate import annotate
from palate_multiome.config import load_config, stage_map
from palate_multiome.io import concat_rna, download_matrix, read_rna
from palate_multiome.qc import apply_gate, compute_metrics
from palate_multiome.rna import embed, normalize, score_cell_cycle


def main(config: str | None) -> None:
    cfg = load_config(config)
    stages = stage_map(cfg)

    adatas, metrics = [], []
    for sample, stage in stages.items():
        path = download_matrix(sample, cfg)
        a = read_rna(path, sample)
        a.obs["stage"] = stage
        adatas.append(a)
        metrics.append(compute_metrics(path, sample, stage, cfg["qc"]["mt_prefix"]))

    adata = concat_rna(adatas)
    del adatas
    adata = apply_gate(adata, pd.concat(metrics, ignore_index=True), cfg["qc"])
    adata.obs["stage"] = pd.Categorical(adata.obs["stage"],
                                        categories=cfg["stage_order"], ordered=True)

    adata = normalize(adata, cfg)
    adata = score_cell_cycle(adata)
    adata = embed(adata, cfg)
    adata = annotate(adata, cfg)

    proc = Path(cfg["data"]["processed_dir"])
    proc.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(proc / "palate_rna.h5ad", compression="gzip")

    out = Path("results")
    out.mkdir(exist_ok=True)
    table = pd.crosstab(adata.obs["celltype"], adata.obs["stage"])
    table.loc["TOTAL"] = table.sum()
    table.to_csv(out / "celltype_by_stage.csv")
    print(table.to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    main(ap.parse_args().config)
