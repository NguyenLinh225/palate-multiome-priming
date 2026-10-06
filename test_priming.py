import numpy as np

from palate_multiome.priming import (
    confident_fraction,
    fate_prediction,
    remove_axis,
    stage_axis,
)


def _separable(n=400, seed=0):
    rng = np.random.default_rng(seed)
    y = np.repeat([0, 1], n // 2)
    X = rng.normal(size=(n, 5))
    X[:, 0] += 3 * y            # fate axis
    return X, y


def test_classifier_recovers_a_planted_signal(cfg):
    X, y = _separable()
    train = np.ones(len(y), bool)
    res = fate_prediction(X, train, y, train, cfg)
    assert res["auroc_mean"] > 0.95
    assert res["frac_confident"] > res["frac_confident_null"]


def test_permuted_labels_give_no_confident_calls(cfg):
    X, y = _separable()
    train = np.ones(len(y), bool)
    res = fate_prediction(X, train, y, train, cfg)
    assert res["frac_confident_null"] < 0.10


def test_pure_noise_yields_no_signal(cfg):
    rng = np.random.default_rng(1)
    X = rng.normal(size=(400, 5))
    y = rng.integers(0, 2, 400)
    train = np.ones(400, bool)
    res = fate_prediction(X, train, y, train, cfg)
    assert res["auroc_mean"] < 0.65


def test_removing_an_axis_destroys_signal_confined_to_it():
    X, y = _separable()
    is_late = y.astype(bool)
    axis = stage_axis(X, is_late, np.ones(len(y), bool))
    Xr = remove_axis(X, axis)
    # the planted separation lived entirely on that axis
    assert abs(Xr[y == 1, 0].mean() - Xr[y == 0, 0].mean()) < 1.0


def test_remove_axis_is_idempotent():
    X, _ = _separable()
    axis = np.zeros(5); axis[0] = 1.0
    once = remove_axis(X, axis)
    twice = remove_axis(once, axis)
    assert np.allclose(once, twice, atol=1e-10)


def test_confident_fraction_respects_the_ambiguous_band(cfg):
    p = np.array([0.01, 0.5, 0.99, 0.6])
    assert confident_fraction(p, cfg) == 0.5
