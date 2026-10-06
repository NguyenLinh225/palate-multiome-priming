"""End-to-end smoke test: run scripts/02_build_rna.py on synthetic 10x files.

Unit tests cover each module; this checks they compose. It runs the real script as a
subprocess against a temporary config, so a broken import, a renamed obs column, or a
shape mismatch between modules fails here even though every unit test passes.
"""
import os
import subprocess
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pytest
import scipy.sparse as sp
import yaml

from palate_multiome.config import load_config

REPO = Path(__file__).resolve().parents[1]
MARKERS = ["Shox2", "Msx1", "Meox2", "Tbx22", "Prrx1", "Twist1", "Col1a1", "Pdgfra",
           "Epcam", "Cdh1", "Krt14", "Runx2", "Sp7", "Sox9", "Acan", "Col2a1"]
# enough of each phase panel that the real cell-cycle scoring path runs
CYCLE = ["Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4",
         "Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Top2a"]


def _write_library(path, sample, n_cells, rng):
    genes = MARKERS + CYCLE + [f"Gene{i}" for i in range(250)] + ["mt-Co1", "mt-Nd1"]
    peaks = [f"chr1:{1000 * i}-{1000 * i + 500}" for i in range(60)]
    rna = rng.poisson(0.6, size=(len(genes), n_cells))
    # two populations so clustering and HVG selection have structure to find
    half = n_cells // 2
    rna[:4, :half] += rng.poisson(8, size=(4, half))
    rna[4:8, half:] += rng.poisson(8, size=(4, n_cells - half))
    atac = rng.poisson(6, size=(len(peaks), n_cells))       # ~360 fragments, clears the gate
    X = sp.csc_matrix(np.vstack([rna, atac]).astype(np.int32))
    with h5py.File(path, "w") as f:
        m = f.create_group("matrix")
        m["data"], m["indices"], m["indptr"] = X.data, X.indices, X.indptr
        m["shape"] = np.array(X.shape, dtype=np.int32)
        m["barcodes"] = np.array([f"BC{i:05d}-1".encode() for i in range(n_cells)])
        ft = m.create_group("features")
        ft["name"] = np.array([g.encode() for g in genes + peaks])
        ft["id"] = np.array([f"G{i}".encode() for i in range(len(genes))]
                            + [p.encode() for p in peaks])
        ft["feature_type"] = np.array([b"Gene Expression"] * len(genes) + [b"Peaks"] * len(peaks))


@pytest.fixture
def synthetic_project(tmp_path):
    rng = np.random.default_rng(0)
    cfg = load_config()
    cfg["samples"] = [{"id": "S1", "stage": "E12.5"}, {"id": "S2", "stage": "E14.5"}]
    raw = tmp_path / "data" / "raw"
    raw.mkdir(parents=True)
    for s in cfg["samples"]:
        _write_library(raw / f"{s['id']}.h5", s["id"], 150, rng)
    cfg["data"]["raw_dir"] = str(raw)
    cfg["data"]["processed_dir"] = str(tmp_path / "data" / "processed")
    cfg["data"]["geo_suppl_url"] = "http://invalid.localhost/"   # must never be contacted
    cfg["rna"].update(n_top_genes=120, n_pcs=20, neighbors_n_pcs=15, neighbors_k=10)
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))
    return tmp_path, cfg_path


def test_build_rna_runs_end_to_end(synthetic_project):
    root, cfg_path = synthetic_project
    # Put src/ on the path explicitly so the subprocess finds the package however
    # pytest was launched (installed, editable, or straight from a checkout).
    env = {**os.environ,
           "PYTHONPATH": os.pathsep.join(filter(None, [str(REPO / "src"),
                                                       os.environ.get("PYTHONPATH")]))}
    proc = subprocess.run(
        [sys.executable, str(REPO / "scripts" / "02_build_rna.py"), "--config", str(cfg_path)],
        cwd=root, env=env, capture_output=True, text=True, timeout=300, check=False,
    )
    assert proc.returncode == 0, proc.stderr[-2000:]

    out = root / "data" / "processed" / "palate_rna.h5ad"
    assert out.exists()
    a = ad.read_h5ad(out)

    # every cell passes the gates by construction, and barcodes collide across libraries
    assert a.n_obs == 300
    assert a.obs_names.is_unique
    for key in ["X_pca", "X_umap"]:
        assert key in a.obsm
    for col in ["sample", "stage", "leiden", "celltype", "nCount_ATAC", "phase"]:
        assert col in a.obs
    assert a.obs["leiden"].nunique() >= 2
    assert set(a.obs["phase"]) <= {"G1", "S", "G2M"}   # real scoring path ran
    assert (root / "results" / "celltype_by_stage.csv").exists()


def test_missing_cell_cycle_genes_degrade_gracefully():
    """A panel without cell-cycle genes must not abort the pipeline."""
    import warnings

    from palate_multiome.rna import score_cell_cycle

    a = ad.AnnData(np.ones((4, 3), dtype=np.float32))
    a.var_names = ["Shox2", "Meox2", "Gapdh"]
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        score_cell_cycle(a)
    assert any("skipping cell-cycle scoring" in str(x.message) for x in w)
    assert (a.obs["phase"] == "unscored").all()
