#!/usr/bin/env python
"""Step 3: the RNA baseline for the priming test.

Trains an anterior-vs-posterior classifier on terminal (E14.5) cells and applies
it to E12.5 progenitors that carry no A/P marker signature. This establishes how
much fate information RNA alone already contains; the ATAC comparison must beat
this number, not merely be non-zero.

Writes results/priming_pilot.csv.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc

from palate_multiome.config import load_config
from palate_multiome.priming import fate_prediction, remove_axis, stage_axis, summarize


def main(config: str | None, h5ad: str | None) -> None:
    cfg = load_config(config)
    p = cfg["priming"]
    path = h5ad or Path(cfg["data"]["processed_dir"]) / "palate_rna.h5ad"
    adata = sc.read_h5ad(path)
    obs = adata.obs

    features = adata.obsm["X_pca"]
    train = ((obs.stage == p["train_stage"]) & obs.celltype.isin(p["classes"])).to_numpy()
    labels = (obs.celltype[train] == p["classes"][0]).astype(int).to_numpy()
    apply_to = ((obs.stage == p["apply_stage"])
                & (obs.celltype == p["apply_celltype"])).to_numpy()

    if p["remove_stage_axis"]:
        both = (obs.stage.isin([p["apply_stage"], p["train_stage"]])
                & obs.celltype.isin(p["stage_axis_celltypes"])).to_numpy()
        is_late = (obs.stage == p["train_stage"]).to_numpy()
        features = remove_axis(features, stage_axis(features, is_late, both, cfg["seed"]))

    res = fate_prediction(features, train, labels, apply_to, cfg)
    out = Path("results")
    out.mkdir(exist_ok=True)
    table = summarize(res, "RNA principal components")
    table.insert(0, "train_set",
                 f"{p['train_stage']} {p['classes'][0]} vs {p['classes'][1]} "
                 f"(n={res['n_train']}, AUROC {res['auroc_mean']:.3f})")
    table.insert(1, "apply_to", f"{p['apply_stage']} {p['apply_celltype']}")
    table.to_csv(out / "priming_pilot.csv", index=False)

    # Ablation: the same test without the stage-axis and class-balance controls, so the
    # effect of each control on the permuted-label null is on record.
    raw_features = adata.obsm["X_pca"]
    naive_cfg = {**cfg, "priming": {**p, "class_weight": None}}
    naive = fate_prediction(raw_features, train, labels, apply_to, naive_cfg)
    ablation = pd.DataFrame([
        {"arm": "naive (no stage-axis removal, unbalanced)",
         "frac_confident": round(naive["frac_confident"], 3),
         "frac_confident_null": round(naive["frac_confident_null"], 3)},
        {"arm": "with controls",
         "frac_confident": round(res["frac_confident"], 3),
         "frac_confident_null": round(res["frac_confident_null"], 3)},
    ])
    ablation.to_csv(out / "priming_controls_ablation.csv", index=False)

    print(f"training AUROC {res['auroc_mean']:.3f} +/- {res['auroc_sd']:.3f} "
          f"(n={res['n_train']})")
    print(f"confident calls on {res['n_apply']:,} progenitors: "
          f"{res['frac_confident']:.1%} real vs {res['frac_confident_null']:.1%} permuted")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--h5ad", default=None)
    a = ap.parse_args()
    main(a.config, a.h5ad)
