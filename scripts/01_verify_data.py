#!/usr/bin/env python
"""Step 1: confirm GEO has usable processed files, and measure the peak-set problem.

Run this before anything else. It downloads only the matrices (~510 MB), not the
fragments (12.3 GB), and answers two questions:
  1. do all libraries carry both modalities?
  2. do the per-library peak sets share a feature space? (they do not)

Writes results/day1_verification.csv and results/peak_id_sharing.csv.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from palate_multiome.config import load_config, stage_map
from palate_multiome.io import download_matrix, read_peak_ids
from palate_multiome.peaks import merge_intervals, parse_intervals, sharing_profile
from palate_multiome.qc import compute_metrics, gate_summary


def main(config: str | None) -> None:
    cfg = load_config(config)
    stages = stage_map(cfg)
    out = Path("results")
    out.mkdir(exist_ok=True)

    metrics, peak_sets, peak_counts = [], {}, {}
    for sample, stage in stages.items():
        path = download_matrix(sample, cfg)
        metrics.append(compute_metrics(path, sample, stage, cfg["qc"]["mt_prefix"]))
        ids = read_peak_ids(path)
        peak_sets[sample] = set(ids.tolist())
        peak_counts[sample] = len(ids)
        print(f"{sample}: {peak_counts[sample]:,} peaks")

    metrics = pd.concat(metrics, ignore_index=True)
    metrics.to_csv(out / "qc_metrics_prefilter.csv.gz", index=False)

    summary = gate_summary(metrics, cfg["qc"])
    summary["peaks_called"] = summary["sample"].map(peak_counts)
    summary.to_csv(out / "day1_verification.csv", index=False)

    sharing = sharing_profile(peak_sets)
    sharing.to_csv(out / "peak_id_sharing.csv", index=False)

    union = set().union(*peak_sets.values())
    consensus = merge_intervals(parse_intervals(union))
    consensus.to_csv(out / "consensus_peaks.bed.gz", sep="\t", index=False,
                     columns=["chrom", "start", "end", "peak_id"], header=False)

    total, passing = summary.input_cells.sum(), summary.pass_count_qc.sum()
    recurring = int(sharing.loc[sharing.n_libraries_sharing_id > 1, "n_peak_ids"].sum())
    print(f"\ncells: {passing:,} of {total:,} pass count QC ({100*passing/total:.1f}%)")
    print(f"peak IDs: {len(union):,} distinct; {recurring:,} recur across libraries "
          f"({100*recurring/len(union):.3f}%)")
    print(f"consensus intervals after merge: {len(consensus):,}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    main(ap.parse_args().config)
