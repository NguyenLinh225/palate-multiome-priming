import numpy as np
import pandas as pd

from palate_multiome.qc import count_gate, gate_summary


def _metrics(**overrides):
    base = {"cell_id": ["c1"], "sample": ["s1"], "stage": ["E12.5"],
            "nCount_RNA": [5000], "nFeature_RNA": [2000], "percent_mt": [5.0],
            "nCount_ATAC": [8000]}
    base.update({k: [v] for k, v in overrides.items()})
    return pd.DataFrame(base)


def test_clean_cell_passes(cfg):
    assert count_gate(_metrics(), cfg["qc"]).iloc[0]


def test_each_gate_rejects_independently(cfg):
    qc = cfg["qc"]
    for field, bad in [("nCount_RNA", qc["max_count_rna"] + 1),
                       ("nFeature_RNA", qc["max_feature_rna"] + 1),
                       ("percent_mt", qc["max_percent_mt"] + 1),
                       ("nCount_ATAC", qc["min_count_atac"] - 1)]:
        assert not count_gate(_metrics(**{field: bad}), qc).iloc[0], field


def test_atac_upper_bound_rejects_doublet_like_cells(cfg):
    bad = _metrics(nCount_ATAC=cfg["qc"]["max_count_atac"] + 1)
    assert not count_gate(bad, cfg["qc"]).iloc[0]


def test_gate_summary_counts_and_percentages(cfg):
    m = pd.concat([_metrics(cell_id=f"c{i}") for i in range(3)]
                  + [_metrics(cell_id="bad", percent_mt=99.0)], ignore_index=True)
    out = gate_summary(m, cfg["qc"])
    assert out.loc[0, "input_cells"] == 4
    assert out.loc[0, "pass_count_qc"] == 3
    assert np.isclose(out.loc[0, "pct_pass"], 75.0)
