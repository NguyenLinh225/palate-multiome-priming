"""End-to-end pipeline. Run with: snakemake -c4"""

configfile: "config/config.yaml"

SAMPLES = [s["id"] for s in config["samples"]]
PROCESSED = config["data"]["processed_dir"]


rule all:
    input:
        "results/day1_verification.csv",
        "results/peak_id_sharing.csv",
        "results/celltype_by_stage.csv",
        "results/priming_pilot.csv",
        "figures/fig1_day1_feasibility.png",
        "figures/fig2_rna_umap.png",
        "figures/fig3_story.png",
        "docs/status_deck.pptx",


rule verify_data:
    """Download matrices, reproduce published QC, measure peak-set fragmentation."""
    output:
        "results/day1_verification.csv",
        "results/peak_id_sharing.csv",
        "results/qc_metrics_prefilter.csv.gz",
        "results/consensus_peaks.bed.gz",
    log:
        "logs/01_verify_data.log",
    shell:
        "python scripts/01_verify_data.py --config config/config.yaml > {log} 2>&1"


rule build_rna:
    """QC, normalize, embed, cluster and annotate the RNA modality."""
    input:
        "results/qc_metrics_prefilter.csv.gz",
    output:
        f"{PROCESSED}/palate_rna.h5ad",
        "results/celltype_by_stage.csv",
    log:
        "logs/02_build_rna.log",
    shell:
        "python scripts/02_build_rna.py --config config/config.yaml > {log} 2>&1"


rule priming_rna:
    """RNA baseline for the fate-prediction test."""
    input:
        f"{PROCESSED}/palate_rna.h5ad",
    output:
        "results/priming_pilot.csv",
    log:
        "logs/03_priming_rna.log",
    shell:
        "python scripts/03_priming_rna.py --config config/config.yaml > {log} 2>&1"


rule figures:
    input:
        "results/day1_verification.csv",
        f"{PROCESSED}/palate_rna.h5ad",
        "results/priming_pilot.csv",
    output:
        "figures/fig1_day1_feasibility.png",
        "figures/fig2_rna_umap.png",
        "figures/fig3_story.png",
    log:
        "logs/04_figures.log",
    shell:
        "python scripts/04_figures.py --config config/config.yaml > {log} 2>&1"


rule deck:
    """7-slide project summary built from the regenerated figures."""
    input:
        "figures/fig1_day1_feasibility.png",
        "figures/fig2_rna_umap.png",
        "figures/fig3_story.png",
    output:
        "docs/status_deck.pptx",
    log:
        "logs/05_build_deck.log",
    shell:
        "python scripts/05_build_deck.py > {log} 2>&1"


# Not yet implemented -- the ATAC half of the project.
# rule consensus_matrix:
#     """Re-quantify all cells against consensus_peaks.bed.gz (12.3 GB fragments)."""
# rule priming_atac:
#     """Same test as priming_rna, on LSI components. This is the actual comparison."""
