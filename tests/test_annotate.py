import pandas as pd

from palate_multiome.annotate import call_lineage, call_mesenchymal_subtype


def test_mesenchymal_cluster_is_not_called_epithelial(cfg):
    scores = pd.DataFrame({"Epithelial": [0.9], "Mesenchymal": [1.8]}, index=["0"])
    assert call_lineage(scores, cfg["lineage_thresholds"])["0"] == "Mesenchymal"


def test_specific_lineage_wins_over_mesenchymal(cfg):
    scores = pd.DataFrame({"Erythroid": [4.2], "Mesenchymal": [1.3]}, index=["0"])
    assert call_lineage(scores, cfg["lineage_thresholds"])["0"] == "Erythroid"


def test_weak_cluster_falls_through_to_other(cfg):
    scores = pd.DataFrame({"Epithelial": [0.1], "Mesenchymal": [0.2]}, index=["0"])
    assert call_lineage(scores, cfg["lineage_thresholds"])["0"] == "Other"


def test_axis_call_requires_dominance_not_just_presence(cfg):
    # anterior above threshold but posterior nearly as high -> unresolved
    scores = pd.DataFrame({"Anterior": [0.9], "Posterior": [0.8], "Osteogenic": [0.1]},
                          index=["0"])
    out = call_mesenchymal_subtype(scores, cfg["subtype_thresholds"])
    assert out["0"] == "Early/unpatterned mesenchyme"


def test_clear_anterior_and_posterior_are_called(cfg):
    scores = pd.DataFrame({"Anterior": [1.5, 0.05], "Posterior": [0.05, 1.3],
                           "Osteogenic": [0.1, 0.1]}, index=["0", "1"])
    out = call_mesenchymal_subtype(scores, cfg["subtype_thresholds"])
    assert out["0"] == "Anterior" and out["1"] == "Posterior"


def test_chondrogenic_requires_acan_detection(cfg):
    """Sox9/Col2a1 without Acan is uncommitted mesenchyme, not cartilage."""
    scores = pd.DataFrame({"Anterior": [0.1], "Posterior": [0.1], "Osteogenic": [0.1],
                           "Chondrogenic": [0.6]}, index=["0"])
    no_acan = call_mesenchymal_subtype(scores, cfg["subtype_thresholds"],
                                       pd.Series({"0": 0.01}))
    with_acan = call_mesenchymal_subtype(scores, cfg["subtype_thresholds"],
                                         pd.Series({"0": 0.55}))
    assert no_acan["0"] == "Early/unpatterned mesenchyme"
    assert with_acan["0"] == "Chondrogenic"
