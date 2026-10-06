"""Reading 10x multiome HDF5 matrices from GEO."""
from __future__ import annotations

import time
import urllib.error
import urllib.request
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

MATRIX_SUFFIX = "_filtered_feature_bc_matrix.h5"
FRAGMENT_SUFFIX = "_atac_fragments.tsv.gz"


def matrix_url(sample: str, cfg: dict) -> str:
    return f"{cfg['data']['geo_suppl_url']}{cfg['data']['accession']}_{sample}{MATRIX_SUFFIX}"


def download_matrix(sample: str, cfg: dict, overwrite: bool = False,
                    retries: int = 3, backoff: float = 5.0) -> Path:
    """Fetch one library's filtered feature-barcode matrix into the raw data dir.

    A file already present at ``data/raw/<sample>.h5`` is used as-is, so the
    matrices can be downloaded manually (browser, ``wget``, Aspera) and dropped in
    place when the NCBI FTP host is unreachable from the current network.
    """
    out_dir = Path(cfg["data"]["raw_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"{sample}.h5"
    if dest.exists() and not overwrite:
        return dest

    url = matrix_url(sample, cfg)
    tmp = dest.with_suffix(".h5.part")
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            urllib.request.urlretrieve(url, tmp)
            tmp.replace(dest)
            return dest
        except urllib.error.URLError as err:
            last_err = err
            if attempt < retries:
                time.sleep(backoff * attempt)
    tmp.unlink(missing_ok=True)
    raise RuntimeError(
        f"could not download {url} after {retries} attempts ({last_err}).\n"
        f"Download it manually and save it as {dest} -- the pipeline will pick it up."
    ) from last_err


def fragment_url(sample: str, cfg: dict) -> str:
    return f"{cfg['data']['geo_suppl_url']}{cfg['data']['accession']}_{sample}{FRAGMENT_SUFFIX}"


def _read_h5(path: str | Path):
    with h5py.File(path) as f:
        m = f["matrix"]
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]),
                          shape=tuple(m["shape"][:]))
        feature_type = m["features"]["feature_type"][:]
        names = np.array([x.decode() for x in m["features"]["name"][:]])
        ids = np.array([x.decode() for x in m["features"]["id"][:]])
        barcodes = np.array([x.decode() for x in m["barcodes"][:]])
    return X, feature_type, names, ids, barcodes


def cell_id(barcode: str, sample: str) -> str:
    """Barcodes collide across libraries, so every cell is suffixed with its sample."""
    return f"{barcode.replace('-1', '')}_{sample}"


def read_rna(path: str | Path, sample: str) -> ad.AnnData:
    """Load the Gene Expression block of a 10x multiome matrix as cells x genes."""
    X, feature_type, names, ids, barcodes = _read_h5(path)
    g = feature_type == b"Gene Expression"
    obs = pd.DataFrame(index=[cell_id(b, sample) for b in barcodes])
    var = pd.DataFrame({"gene_ids": ids[g]}, index=pd.Index(names[g]))
    adata = ad.AnnData(X[g].T.tocsr().astype(np.float32), obs=obs, var=var)
    adata.obs["sample"] = sample
    return adata


def read_peak_ids(path: str | Path) -> np.ndarray:
    """Peak interval IDs (``chr:start-end``) for one library."""
    _, feature_type, _, ids, _ = _read_h5(path)
    return ids[feature_type == b"Peaks"]


def concat_rna(adatas: list[ad.AnnData]) -> ad.AnnData:
    """Stack libraries sharing an identical gene order.

    cellranger-arc emits the same reference gene list for every library, so the
    matrices are stacked directly. The gene order is asserted rather than assumed;
    a mismatch means the libraries were built against different references and
    must be reindexed instead.
    """
    ref = adatas[0].var["gene_ids"].to_numpy()
    for a in adatas[1:]:
        if not np.array_equal(a.var["gene_ids"].to_numpy(), ref):
            raise ValueError("gene order differs between libraries; reindex before concat")
    X = sp.vstack([a.X for a in adatas], format="csr")
    obs = pd.concat([a.obs for a in adatas])
    out = ad.AnnData(X, obs=obs, var=adatas[0].var.copy())
    out.var_names_make_unique()
    return out
