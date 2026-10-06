import h5py
import numpy as np
import pytest
import scipy.sparse as sp

from palate_multiome.io import cell_id, concat_rna, download_matrix, read_peak_ids, read_rna


def _write_10x_h5(path, genes, peaks, barcodes, dense):
    """Minimal cellranger-arc layout: features x barcodes, CSC."""
    X = sp.csc_matrix(np.asarray(dense, dtype=np.int32))
    names = genes + peaks
    ftype = [b"Gene Expression"] * len(genes) + [b"Peaks"] * len(peaks)
    ids = [f"ENSMUSG{i:05d}".encode() for i in range(len(genes))] + [p.encode() for p in peaks]
    with h5py.File(path, "w") as f:
        m = f.create_group("matrix")
        m["data"], m["indices"], m["indptr"] = X.data, X.indices, X.indptr
        m["shape"] = np.array(X.shape, dtype=np.int32)
        m["barcodes"] = np.array([b.encode() for b in barcodes])
        ft = m.create_group("features")
        ft["name"] = np.array([n.encode() for n in names])
        ft["id"] = np.array(ids)
        ft["feature_type"] = np.array(ftype)


@pytest.fixture
def lib(tmp_path):
    def make(name, barcodes=("AAAC-1", "GGTT-1")):
        p = tmp_path / f"{name}.h5"
        # 2 genes + 1 peak, by 2 cells
        _write_10x_h5(p, ["Shox2", "Meox2"], ["chr1:100-200"], list(barcodes),
                      [[5, 0], [0, 3], [2, 1]])
        return p
    return make


def test_cell_id_suffixes_sample_and_strips_gem_suffix():
    assert cell_id("AAAC-1", "B6E12-5-1") == "AAAC_B6E12-5-1"


def test_read_rna_keeps_only_gene_expression_as_cells_by_genes(lib):
    a = read_rna(lib("s1"), "s1")
    assert a.shape == (2, 2)
    assert list(a.var_names) == ["Shox2", "Meox2"]
    assert a.X[0, 0] == 5 and a.X[1, 1] == 3


def test_read_peak_ids_returns_only_peaks(lib):
    assert read_peak_ids(lib("s1")).tolist() == ["chr1:100-200"]


def test_identical_barcodes_in_two_libraries_stay_distinct(lib):
    """10x barcodes recur across libraries; the sample suffix must disambiguate."""
    a = concat_rna([read_rna(lib("s1"), "s1"), read_rna(lib("s2"), "s2")])
    assert a.n_obs == 4
    assert a.obs_names.is_unique


def test_concat_rejects_mismatched_gene_order(lib, tmp_path):
    p = tmp_path / "bad.h5"
    _write_10x_h5(p, ["Meox2", "Shox2"], ["chr1:100-200"], ["AAAC-1", "GGTT-1"],
                  [[5, 0], [0, 3], [2, 1]])
    a = read_rna(lib("s1"), "s1")
    b = read_rna(p, "bad")
    b.var["gene_ids"] = b.var["gene_ids"].iloc[::-1].to_numpy()
    with pytest.raises(ValueError, match="gene order"):
        concat_rna([a, b])


def test_existing_file_is_used_without_network(tmp_path, cfg, lib):
    """Manually downloaded matrices must short-circuit the download."""
    local = dict(cfg)
    local["data"] = {**cfg["data"], "raw_dir": str(tmp_path / "raw"),
                     "geo_suppl_url": "http://invalid.localhost/"}
    (tmp_path / "raw").mkdir()
    lib("s1").rename(tmp_path / "raw" / "s1.h5")
    assert download_matrix("s1", local) == tmp_path / "raw" / "s1.h5"


def test_unreachable_host_gives_actionable_error(tmp_path, cfg):
    local = dict(cfg)
    local["data"] = {**cfg["data"], "raw_dir": str(tmp_path / "raw"),
                     "geo_suppl_url": "http://invalid.localhost/"}
    with pytest.raises(RuntimeError, match="Download it manually"):
        download_matrix("s1", local, retries=1, backoff=0)
    assert not (tmp_path / "raw" / "s1.h5.part").exists()
