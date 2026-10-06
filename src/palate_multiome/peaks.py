"""Consensus peak construction.

The nine libraries were peak-called independently: 67k-140k peaks each, and only
0.08% of peak IDs recur in any other library (none in all nine). The matrices
therefore share no ATAC feature space and cannot be concatenated as deposited.
This module merges the per-library intervals into one consensus set that all
cells can be re-quantified against.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd


def parse_intervals(peak_ids) -> list[tuple[str, int, int]]:
    """``chr1:3669052-3670290`` -> ``('chr1', 3669052, 3670290)``."""
    out = []
    for pid in peak_ids:
        pid = pid.decode() if isinstance(pid, bytes) else pid
        chrom, _, rest = pid.partition(":")
        start, _, end = rest.partition("-")
        out.append((chrom, int(start), int(end)))
    return out


def merge_intervals(intervals) -> pd.DataFrame:
    """Union-merge overlapping intervals into a sorted consensus set."""
    by_chrom = defaultdict(list)
    for chrom, start, end in intervals:
        by_chrom[chrom].append((start, end))

    rows = []
    for chrom in sorted(by_chrom):
        spans = sorted(by_chrom[chrom])
        cur_start, cur_end = spans[0]
        for start, end in spans[1:]:
            if start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                rows.append((chrom, cur_start, cur_end))
                cur_start, cur_end = start, end
        rows.append((chrom, cur_start, cur_end))

    df = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    df["width"] = df["end"] - df["start"]
    df["peak_id"] = df.chrom + ":" + df.start.astype(str) + "-" + df.end.astype(str)
    return df


def sharing_profile(peak_sets: dict[str, set]) -> pd.DataFrame:
    """How many peak IDs are shared by exactly N libraries.

    The motivating diagnostic: in GSE218576 this returns 895,069 IDs in one
    library, 727 in two, 1 in three, from a 895,797-ID union.
    """
    counts: dict[str, int] = defaultdict(int)
    for ids in peak_sets.values():
        for pid in ids:
            counts[pid] += 1
    tally = pd.Series(list(counts.values())).value_counts().sort_index()
    return pd.DataFrame({"n_libraries_sharing_id": tally.index, "n_peak_ids": tally.to_numpy()})


def filter_by_prevalence(counts, min_cell_fraction: float) -> np.ndarray:
    """Boolean mask of peaks detected in at least ``min_cell_fraction`` of cells."""
    n_cells = counts.shape[0]
    detected = np.asarray((counts > 0).sum(0)).ravel()
    return detected >= min_cell_fraction * n_cells
