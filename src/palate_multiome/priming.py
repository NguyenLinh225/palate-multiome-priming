"""The fate-prediction test.

Train a classifier to separate terminal anterior from posterior cells, then
apply it to progenitors that carry no anterior/posterior marker signature. The
question is whether an earlier cell's state already predicts where it ends up.

Running this on RNA features gives the baseline the chromatin claim must beat;
the same function on ATAC features gives the comparison. Three controls are
built in, because each corresponds to a way the naive result can be spurious:

* ``remove_stage_axis`` projects out the E12.5-vs-E14.5 discriminant, so
  apparent fate signal cannot be maturation signal.
* ``class_weight='balanced'`` prevents the majority class from producing
  confident calls by base rate alone.
* ``permutation_null`` reruns with shuffled training labels.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score


def stage_axis(features: np.ndarray, is_late: np.ndarray, subset: np.ndarray,
               seed: int = 0) -> np.ndarray:
    """Unit vector separating early from late cells, to be projected out."""
    clf = LogisticRegression(max_iter=2000, random_state=seed)
    clf.fit(features[subset], is_late[subset])
    w = clf.coef_[0]
    return w / np.linalg.norm(w)


def remove_axis(features: np.ndarray, axis: np.ndarray) -> np.ndarray:
    """Project features orthogonal to ``axis``."""
    return features - np.outer(features @ axis, axis)


def fate_prediction(features: np.ndarray, train_mask: np.ndarray, train_labels: np.ndarray,
                    apply_mask: np.ndarray, cfg: dict, seed: int | None = None) -> dict:
    """Train on terminal cells, score progenitors, and compare to a shuffled null.

    Parameters
    ----------
    features
        Cells x components. RNA PCs for the baseline, LSI components for ATAC.
    train_labels
        1 for the positive class (anterior), 0 for the negative.

    Returns
    -------
    dict with cross-validated training AUROC, the progenitor probabilities, the
    fraction receiving a confident call, and the same under permuted labels.
    """
    p = cfg["priming"]
    seed = cfg["seed"] if seed is None else seed
    rng = np.random.default_rng(seed)

    def _fit(y):
        return LogisticRegression(max_iter=2000, class_weight=p["class_weight"],
                                  random_state=seed).fit(features[train_mask], y)

    cv = cross_val_score(
        LogisticRegression(max_iter=2000, class_weight=p["class_weight"], random_state=seed),
        features[train_mask], train_labels,
        cv=StratifiedKFold(p["cv_folds"], shuffle=True, random_state=seed),
        scoring="roc_auc")

    proba = _fit(train_labels).predict_proba(features[apply_mask])[:, 1]
    null = _fit(rng.permutation(train_labels)).predict_proba(features[apply_mask])[:, 1]

    return {
        "auroc_mean": float(cv.mean()),
        "auroc_sd": float(cv.std()),
        "n_train": int(train_mask.sum()),
        "n_apply": int(apply_mask.sum()),
        "proba": proba,
        "proba_null": null,
        "frac_confident": confident_fraction(proba, cfg),
        "frac_confident_null": confident_fraction(null, cfg),
    }


def confident_fraction(proba: np.ndarray, cfg: dict) -> float:
    """Fraction of cells assigned a fate with probability outside the ambiguous band."""
    p = cfg["priming"]
    return float(np.mean((proba < p["confident_low"]) | (proba > p["confident_high"])))


def summarize(result: dict, label: str) -> pd.DataFrame:
    """One row per arm, for writing to results/."""
    rows = []
    for arm, key in [("real labels", "proba"), ("permuted labels", "proba_null")]:
        v = result[key]
        rows.append({
            "features": label,
            "arm": arm,
            "n": len(v),
            "mean_p": round(float(v.mean()), 3),
            "sd": round(float(v.std()), 3),
            "frac_confident": round(float(np.mean((v < 0.2) | (v > 0.8))), 3),
        })
    return pd.DataFrame(rows)
