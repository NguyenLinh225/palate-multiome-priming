#!/usr/bin/env python
"""Step 4: regenerate the three figures from the committed results and the RNA object."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.lines import Line2D

from palate_multiome.config import load_config

CELLTYPE_COLORS = {
    "Anterior": "#1B7837", "Posterior": "#762A83", "Osteogenic": "#D9822B",
    "Early/unpatterned mesenchyme": "#C7C7C7", "Epithelial": "#2166AC",
    "Neural/glia": "#8C6BB1", "Myocyte": "#66A61E", "Endothelial": "#E7298A",
    "Immune": "#A6761D", "Erythroid": "#E31A1C", "Other": "#E8E8E8",
}
RUST = "#B3541E"
GREY = "#7A7A7A"


def style() -> None:
    mpl.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 300, "font.size": 8,
        "axes.titlesize": 8, "axes.labelsize": 8, "axes.spines.top": False,
        "axes.spines.right": False, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "legend.fontsize": 6, "axes.titlelocation": "left",
    })


def panel_letter(ax, letter: str) -> None:
    ax.annotate(letter, xy=(0, 1), xytext=(-24, 10), xycoords="axes fraction",
                textcoords="offset points", fontsize=11, fontweight="bold", va="top")


def bare(ax) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def figure1(cfg, qc, verification, cpm, out: Path) -> None:
    """QC distributions, retention per library, peak-set fragmentation, markers."""
    stages = cfg["stage_order"]
    pal = dict(zip(stages, plt.cm.BuGn([0.42, 0.58, 0.72, 0.9])))
    fig = plt.figure(figsize=(7.2, 6.8))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.05], hspace=0.68, wspace=0.40)

    specs = [("nFeature_RNA", "Genes per cell", cfg["qc"]["max_feature_rna"], None),
             ("percent_mt", "Mitochondrial reads (%)", cfg["qc"]["max_percent_mt"], None),
             ("nCount_ATAC", "ATAC fragments per cell", cfg["qc"]["min_count_atac"], "log")]
    first = None
    for j, (col, lab, thr, scale) in enumerate(specs):
        ax = fig.add_subplot(gs[0, j])
        first = first or ax
        data = [qc.loc[qc.stage == s, col].to_numpy() for s in stages]
        vp = ax.violinplot(data, showextrema=False, widths=0.85)
        for body, s in zip(vp["bodies"], stages):
            body.set_facecolor(pal[s])
            body.set_alpha(0.9)
            body.set_edgecolor("none")
        ax.boxplot(data, widths=0.13, showfliers=False,
                   medianprops=dict(color="white", lw=1.1),
                   boxprops=dict(color=GREY, lw=0.8),
                   whiskerprops=dict(color=GREY, lw=0.8), capprops=dict(lw=0))
        if scale == "log":
            ax.set_yscale("log")
            ax.axhline(cfg["qc"]["max_count_atac"], ls="--", lw=0.9, color=RUST)
        ax.axhline(thr, ls="--", lw=0.9, color=RUST)
        ax.set_xticks(range(1, len(stages) + 1))
        ax.set_xticklabels(stages, rotation=45, ha="right")
        ax.set_ylabel(lab)

    axb = fig.add_subplot(gs[1, :2])
    x = np.arange(len(verification))
    axb.bar(x, verification.input_cells, color="#D9D9D9", width=0.72,
            label="All barcodes in matrix")
    axb.bar(x, verification.pass_count_qc,
            color=[pal[s] for s in verification.stage], width=0.72,
            label="Pass RNA/ATAC-count QC")
    for xi, (inp, passed) in enumerate(zip(verification.input_cells,
                                           verification.pass_count_qc)):
        axb.text(xi, passed * 0.45, f"{100 * passed / inp:.0f}%", ha="center",
                 va="center", fontsize=6, color="white")
    axb.set_xticks(x)
    axb.set_xticklabels([s.replace("B6", "") for s in verification["sample"]],
                        rotation=45, ha="right")
    axb.set_ylabel("Cells")
    axb.legend(frameon=False, loc="upper right")
    total, passed = verification.input_cells.sum(), verification.pass_count_qc.sum()
    axb.set_title(f"{passed:,} of {total:,} barcodes clear the count-based gates")

    axc = fig.add_subplot(gs[1, 2])
    axc.barh(np.arange(len(verification)), verification.peaks_called,
             color=[pal[s] for s in verification.stage], height=0.72)
    axc.set_yticks(np.arange(len(verification)))
    axc.set_yticklabels([s.replace("B6E", "E") for s in verification["sample"]], fontsize=6)
    axc.invert_yaxis()
    axc.set_xlabel("Peaks called")
    axc.set_title("Peak sets are library-specific:\n0.08% of IDs recur, none in all 9",
                  fontsize=7)

    axd = fig.add_subplot(gs[2, :])
    groups = {"Anterior": ["Shox2", "Msx1"], "Posterior": ["Meox2", "Tbx22"],
              "Osteogenic": ["Runx2", "Sp7"]}
    colors = {"Anterior": "#1B7837", "Posterior": "#762A83", "Osteogenic": "#D9822B"}
    xs = np.arange(len(stages))
    ends = {}
    for group, genes in groups.items():
        for k, gene in enumerate(genes):
            y = cpm.loc[stages, gene].to_numpy()
            axd.plot(xs, y, ["-", "--"][k], color=colors[group], lw=1.4, marker="o", ms=3.5)
            ends[gene] = (y[-1], colors[group])
    order = sorted(ends, key=lambda g: -ends[g][0])
    logy = {g: np.log10(max(ends[g][0], 1e-2)) for g in order}
    for i in range(1, len(order)):
        if logy[order[i - 1]] - logy[order[i]] < 0.13:
            logy[order[i]] = logy[order[i - 1]] - 0.13
    for gene in order:
        axd.plot([xs[-1] + 0.015, xs[-1] + 0.085], [ends[gene][0], 10 ** logy[gene]],
                 lw=0.6, color=ends[gene][1], clip_on=False)
        axd.text(xs[-1] + 0.10, 10 ** logy[gene], f"$\\it{{{gene}}}$", fontsize=6.5,
                 color=ends[gene][1], va="center")
    axd.set_xticks(xs)
    axd.set_xticklabels(stages)
    axd.set_xlim(-0.15, len(stages) - 1 + 0.40)
    axd.set_yscale("log")
    axd.set_ylabel("Pseudobulk CPM")
    axd.set_yticks([1, 10, 100, 1000])
    axd.set_yticklabels(["1", "10", "100", "1k"])
    for group, color in colors.items():
        axd.plot([], [], color=color, lw=1.4, label=group)
    axd.legend(frameon=False, loc="upper left", ncol=3, bbox_to_anchor=(0.14, 1.02))
    axd.set_title("Anterior, posterior and osteogenic markers show expected stage dynamics")

    for ax, letter in zip([first, axb, axc, axd], "abcd"):
        panel_letter(ax, letter)
    fig.savefig(out / "fig1_day1_feasibility.png", bbox_inches="tight")
    plt.close(fig)


def figure2(adata, cfg, out: Path) -> None:
    """UMAP coloured by cell type, stage, and the two axis markers."""
    stages = cfg["stage_order"]
    st_cols = dict(zip(stages, plt.cm.BuGn([0.42, 0.58, 0.72, 0.9])))
    U = adata.obsm["X_umap"]
    obs = adata.obs
    source = adata.raw.to_adata()

    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 2, hspace=0.20, wspace=0.08)

    axa = fig.add_subplot(gs[0, 0])
    colors = np.array([CELLTYPE_COLORS.get(v, "#CCCCCC") for v in obs.celltype])
    order = np.argsort([0 if v in ("Early/unpatterned mesenchyme", "Other") else 1
                        for v in obs.celltype])
    axa.scatter(U[order, 0], U[order, 1], c=colors[order], s=1.1, lw=0, rasterized=True)
    bare(axa)
    axa.margins(0.12)
    for ct in ["Anterior", "Posterior", "Osteogenic", "Epithelial",
               "Early/unpatterned mesenchyme"]:
        m = (obs.celltype == ct).to_numpy()
        if not m.any():
            continue
        text = "Early/unpatterned\nmesenchyme" if ct.startswith("Early") else ct
        axa.text(np.median(U[m, 0]), np.median(U[m, 1]), text, fontsize=6.2,
                 ha="center", va="center",
                 bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.0))
    minor = ["Neural/glia", "Myocyte", "Endothelial", "Immune", "Erythroid"]
    axa.legend(handles=[Line2D([], [], marker="o", ls="", ms=4,
                               color=CELLTYPE_COLORS[c], label=c)
                        for c in minor if (obs.celltype == c).any()],
               frameon=False, loc="lower right")
    axa.set_title("Transcriptional identity")

    axb = fig.add_subplot(gs[0, 1])
    perm = np.random.default_rng(cfg["seed"]).permutation(len(U))
    axb.scatter(U[perm, 0], U[perm, 1], c=[st_cols[v] for v in obs.stage[perm]],
                s=1.1, lw=0, rasterized=True)
    bare(axb)
    axb.margins(0.12)
    axb.legend(handles=[Line2D([], [], marker="o", ls="", ms=4, color=st_cols[s], label=s)
                        for s in stages], frameon=False, loc="lower right")
    axb.set_title("Developmental stage")

    axes = []
    for j, (gene, color, side) in enumerate([("Shox2", "#1B7837", "anterior"),
                                             ("Meox2", "#762A83", "posterior")]):
        ax = fig.add_subplot(gs[1, j])
        axes.append(ax)
        v = np.asarray(source[:, [gene]].X.todense()).ravel()
        o = np.argsort(v)
        cmap = mpl.colors.LinearSegmentedColormap.from_list(gene, ["#EAEAEA", color])
        h = ax.scatter(U[o, 0], U[o, 1], c=v[o], s=1.1, lw=0, cmap=cmap,
                       rasterized=True, vmin=0, vmax=np.quantile(v, 0.995))
        bare(ax)
        ax.margins(0.12)
        cb = fig.colorbar(h, ax=ax, fraction=0.028, pad=0.005)
        cb.set_label("log-norm expression", fontsize=6)
        cb.ax.tick_params(labelsize=6)
        cb.outline.set_visible(False)
        ax.set_title(f"$\\it{{{gene}}}$, {side} marker")

    for ax, letter in zip([axa, axb, axes[0], axes[1]], "abcd"):
        panel_letter(ax, letter)
    fig.savefig(out / "fig2_rna_umap.png", bbox_inches="tight")
    plt.close(fig)


def figure3(adata, cfg, proba, proba_null, out: Path) -> None:
    """Patterning onset, the test window, and the fate-prediction result."""
    stages = cfg["stage_order"]
    obs = adata.obs
    U = adata.obsm["X_umap"]
    mes = ["Anterior", "Posterior", "Osteogenic", "Early/unpatterned mesenchyme"]
    sub = obs[obs.celltype.isin(mes)]
    frac = pd.crosstab(sub.stage, sub.celltype, normalize="index") * 100

    fig = plt.figure(figsize=(7.2, 2.6))
    gs = fig.add_gridspec(1, 3, wspace=0.42)

    axa = fig.add_subplot(gs[0, 0])
    bottom = np.zeros(len(stages))
    for ct in ["Early/unpatterned mesenchyme", "Anterior", "Posterior", "Osteogenic"]:
        if ct not in frac:
            continue
        v = frac[ct].reindex(stages).fillna(0).to_numpy()
        axa.bar(np.arange(len(stages)), v, bottom=bottom, color=CELLTYPE_COLORS[ct],
                width=0.72, label="Unpatterned" if ct.startswith("Early") else ct)
        bottom += v
    axa.set_xticks(range(len(stages)))
    axa.set_xticklabels(stages, rotation=45, ha="right")
    axa.set_ylabel("% of mesenchyme")
    axa.set_ylim(0, 100)
    axa.legend(frameon=False, loc="lower left")
    axa.set_title("A/P identity appears at E13.5")

    axb = fig.add_subplot(gs[0, 1])
    early = (obs.stage == cfg["priming"]["apply_stage"]).to_numpy() & obs.celltype.isin(mes).to_numpy()
    late = (obs.stage == cfg["priming"]["train_stage"]).to_numpy() & obs.celltype.isin(mes).to_numpy()
    axb.scatter(U[~(early | late), 0], U[~(early | late), 1], c="#EDEDED", s=0.8,
                lw=0, rasterized=True)
    axb.scatter(U[late, 0], U[late, 1],
                c=[CELLTYPE_COLORS[v] for v in obs.celltype[late]], s=1.0, lw=0,
                rasterized=True)
    axb.scatter(U[early, 0], U[early, 1], c="#3B3B3B", s=1.0, lw=0, rasterized=True)
    bare(axb)
    axb.margins(0.10)
    axb.annotate("E12.5 progenitors\n(no A/P signature)",
                 xy=(np.median(U[early, 0]), np.median(U[early, 1])),
                 xytext=(0.03, 0.06), textcoords="axes fraction", fontsize=6,
                 arrowprops=dict(arrowstyle="-", lw=0.6, color="#3B3B3B"))
    axb.set_title("The test window")

    axc = fig.add_subplot(gs[0, 2])
    bins = np.linspace(0, 1, 36)
    axc.hist(proba_null, bins=bins, color="#BDBDBD", lw=0, label="Shuffled labels")
    axc.hist(proba, bins=bins, color="#1B7837", alpha=0.8, lw=0, label="Real labels")
    axc.axvline(0.5, ls="--", lw=0.8, color=GREY)
    axc.set_xlabel("P(anterior) from RNA")
    axc.set_ylabel(f"{cfg['priming']['apply_stage']} progenitors")
    axc.legend(frameon=False, loc="upper left")
    axc.set_title("RNA already predicts fate")

    for ax, letter in zip([axa, axb, axc], "abc"):
        panel_letter(ax, letter)
    fig.savefig(out / "fig3_story.png", bbox_inches="tight")
    plt.close(fig)


def pseudobulk_cpm(adata, cfg) -> pd.DataFrame:
    """Counts-per-million per stage for the marker genes in panel 1d."""
    genes = ["Shox2", "Msx1", "Meox2", "Tbx22", "Runx2", "Sp7"]
    counts = adata.raw.to_adata() if adata.raw is not None else adata
    rows = {}
    for stage in cfg["stage_order"]:
        m = (adata.obs["stage"] == stage).to_numpy()
        total = np.asarray(counts[m].X.expm1().sum(0)).ravel()
        cpm = total / total.sum() * 1e6
        rows[stage] = {g: cpm[counts.var_names.get_loc(g)]
                       for g in genes if g in counts.var_names}
    return pd.DataFrame(rows).T


def rerun_pilot(adata, cfg):
    """Recompute progenitor fate probabilities so panel 3c can show distributions."""
    from palate_multiome.priming import fate_prediction, remove_axis, stage_axis

    p = cfg["priming"]
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
    return res["proba"], res["proba_null"]


def main(config: str | None) -> None:
    cfg = load_config(config)
    style()
    out = Path("figures")
    out.mkdir(exist_ok=True)
    results = Path("results")

    adata = sc.read_h5ad(Path(cfg["data"]["processed_dir"]) / "palate_rna.h5ad")

    qc = pd.read_csv(results / "qc_metrics_prefilter.csv.gz")
    verification = pd.read_csv(results / "day1_verification.csv")
    figure1(cfg, qc, verification, pseudobulk_cpm(adata, cfg), out)

    figure2(adata, cfg, out)
    proba, proba_null = rerun_pilot(adata, cfg)
    figure3(adata, cfg, proba, proba_null, out)
    print(f"wrote 3 figures to {out}/")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    main(ap.parse_args().config)
