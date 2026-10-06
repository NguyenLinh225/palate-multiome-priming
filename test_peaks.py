import numpy as np
import pandas as pd
import scipy.sparse as sp

from palate_multiome.peaks import (
    filter_by_prevalence,
    merge_intervals,
    parse_intervals,
    sharing_profile,
)


def test_parse_handles_str_and_bytes():
    assert parse_intervals(["chr1:100-200"]) == [("chr1", 100, 200)]
    assert parse_intervals([b"chr2:5-9"]) == [("chr2", 5, 9)]


def test_overlapping_intervals_merge():
    out = merge_intervals([("chr1", 100, 200), ("chr1", 150, 300)])
    assert len(out) == 1
    assert out.loc[0, "start"] == 100 and out.loc[0, "end"] == 300


def test_disjoint_intervals_are_kept_separate():
    out = merge_intervals([("chr1", 100, 200), ("chr1", 500, 600)])
    assert len(out) == 2


def test_chromosomes_do_not_merge_across():
    out = merge_intervals([("chr1", 100, 200), ("chr2", 150, 250)])
    assert len(out) == 2
    assert set(out.chrom) == {"chr1", "chr2"}


def test_merge_is_order_independent():
    a = merge_intervals([("chr1", 100, 200), ("chr1", 150, 300), ("chr1", 900, 950)])
    b = merge_intervals([("chr1", 900, 950), ("chr1", 150, 300), ("chr1", 100, 200)])
    pd.testing.assert_frame_equal(a, b)


def test_sharing_profile_counts_library_recurrence():
    out = sharing_profile({"a": {"p1", "p2"}, "b": {"p2", "p3"}})
    counts = dict(zip(out.n_libraries_sharing_id, out.n_peak_ids))
    assert counts == {1: 2, 2: 1}


def test_prevalence_filter_uses_cell_fraction():
    # peak 0 in 3/4 cells, peak 1 in 1/4
    counts = sp.csr_matrix(np.array([[1, 0], [1, 0], [1, 1], [0, 0]]))
    keep = filter_by_prevalence(counts, 0.5)
    assert keep.tolist() == [True, False]
