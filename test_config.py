import pytest

from palate_multiome.config import load_config, sample_ids, stage_map


def test_loads_and_validates(cfg):
    assert cfg["seed"] == 1234
    assert len(sample_ids(cfg)) == 9


def test_every_sample_has_a_known_stage(cfg):
    assert set(stage_map(cfg).values()) <= set(cfg["stage_order"])


def test_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        load_config("does/not/exist.yaml")


def test_qc_gates_are_ordered(cfg):
    qc = cfg["qc"]
    assert qc["min_count_atac"] < qc["max_count_atac"]
    assert 0 < qc["max_percent_mt"] <= 100
