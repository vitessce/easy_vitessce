import pytest
import scanpy as sc
from spatialdata import SpatialData

import easy_vitessce as ev
from easy_vitessce.data import (
    OBJECT_TO_PARAMS,
    _get_adata_wrapper_params,
    _get_sdata_wrapper_params,
)


@pytest.fixture(autouse=True)
def clear_registry():
    # Registrations are keyed by id(), which can be reused after garbage collection,
    # so clear them to avoid leaking into other tests.
    yield
    OBJECT_TO_PARAMS.clear()


def test_register_data_params_adata_path():
    adata = sc.datasets.pbmc68k_reduced()
    ev.register_data_params(adata, {"adata_path": "/path/to/my_object.adata.zarr"})
    assert _get_adata_wrapper_params(adata) == {"adata_path": "/path/to/my_object.adata.zarr"}


def test_register_data_params_adata_url():
    adata = sc.datasets.pbmc68k_reduced()
    ev.register_data_params(adata, {"adata_url": "https://example.com/my_object.adata.zarr"})
    assert _get_adata_wrapper_params(adata) == {"adata_url": "https://example.com/my_object.adata.zarr"}


def test_register_data_params_adata_store():
    adata = sc.datasets.pbmc68k_reduced()
    store = object()  # Any store-like object should be passed through as-is.
    ev.register_data_params(adata, {"adata_store": store})
    params = _get_adata_wrapper_params(adata)
    assert params == {"adata_store": store}
    assert params["adata_store"] is store


def test_register_data_params_returns_copy():
    adata = sc.datasets.pbmc68k_reduced()
    ev.register_data_params(adata, {"adata_path": "a.adata.zarr"})
    params = _get_adata_wrapper_params(adata)
    params["adata_path"] = "b.adata.zarr"
    assert _get_adata_wrapper_params(adata) == {"adata_path": "a.adata.zarr"}


def test_register_data_params_sdata():
    sdata = SpatialData()
    ev.register_data_params(sdata, {"sdata_url": "https://example.com/my_object.sdata.zarr"})
    assert _get_sdata_wrapper_params(sdata) == {"sdata_url": "https://example.com/my_object.sdata.zarr"}


def test_register_data_params_rejects_non_dict():
    adata = sc.datasets.pbmc68k_reduced()
    with pytest.raises(TypeError):
        ev.register_data_params(adata, "/path/to/my_object.adata.zarr")


def test_register_data_path():
    adata = sc.datasets.pbmc68k_reduced()
    ev.register_data_path(adata, "/path/to/my_object.adata.zarr")
    assert _get_adata_wrapper_params(adata) == {"adata_path": "/path/to/my_object.adata.zarr"}

    sdata = SpatialData()
    ev.register_data_path(sdata, "/path/to/my_object.sdata.zarr")
    assert _get_sdata_wrapper_params(sdata) == {"sdata_path": "/path/to/my_object.sdata.zarr"}


def test_register_data_store():
    adata = sc.datasets.pbmc68k_reduced()
    store = object()
    ev.register_data_store(adata, store)
    assert _get_adata_wrapper_params(adata)["adata_store"] is store

    sdata = SpatialData()
    ev.register_data_store(sdata, store)
    assert _get_sdata_wrapper_params(sdata)["sdata_store"] is store


def test_register_data_url():
    adata = sc.datasets.pbmc68k_reduced()
    ev.register_data_url(adata, "https://example.com/my_object.adata.zarr")
    assert _get_adata_wrapper_params(adata) == {"adata_url": "https://example.com/my_object.adata.zarr"}

    ev.register_data_url(adata, "https://example.com/my_object.h5ad")
    assert _get_adata_wrapper_params(adata) == {
        "adata_url": "https://example.com/my_object.h5ad",
        "ref_url": "https://example.com/my_object.ref.json",
    }

    sdata = SpatialData()
    ev.register_data_url(sdata, "https://example.com/my_object.sdata.zarr")
    assert _get_sdata_wrapper_params(sdata) == {"sdata_url": "https://example.com/my_object.sdata.zarr"}


def test_register_convenience_rejects_unknown_type():
    with pytest.raises(TypeError):
        ev.register_data_path({}, "/path/to/my_object.adata.zarr")


def test_written_adata_is_registered(tmp_path):
    adata = sc.datasets.pbmc68k_reduced()
    ev.config.set({"data.out_dir": str(tmp_path)})
    try:
        params = _get_adata_wrapper_params(adata)
        assert params["adata_path"].startswith(str(tmp_path))
        assert OBJECT_TO_PARAMS[id(adata)] == params
        # Subsequent calls reuse the written file rather than writing a new one.
        assert _get_adata_wrapper_params(adata) == params
        assert len(list(tmp_path.iterdir())) == 1
    finally:
        ev.config.refresh()
