import os
from os.path import join, isfile, isdir, exists
from uuid import uuid4
import shutil
import warnings
import json

from vitessce.data_utils import VAR_CHUNK_SIZE, generate_h5ad_ref_spec

from .widget import config

# AnnData/SpatialData objects are not hashable,
# so we use the id() function to get a unique identifier for the object.
# Reference: https://github.com/scverse/anndata/issues/742

# Maps id(obj) -> dict of data-related wrapper parameters.
# Populated by the user via register_data_params (or its convenience wrappers),
# and by easy_vitessce after writing an object to disk, to avoid writing it again.
OBJECT_TO_PARAMS = {}


def register_data_params(obj, params):
    """
    Register an AnnData or SpatialData object with the data-related parameters
    that should be passed to the corresponding Vitessce data wrapper
    (AnnDataWrapper or SpatialDataWrapper).
    This allows to avoid writing a large object to disk if it has already been written,
    and allows using alternative storage mechanisms (e.g., Zarr stores or remote URLs).

    For AnnData objects, valid keys include ``adata_path``, ``adata_url``, ``adata_store``,
    ``ref_path``, and ``ref_url``.
    For SpatialData objects, valid keys include ``sdata_path``, ``sdata_url``, and ``sdata_store``.

    .. code-block:: python

        ev.register_data_params(adata, { "adata_path": "/path/to/my_object.adata.zarr" })
        ev.register_data_params(adata, { "adata_store": zarr.storage.ZipStore("/path/to/my_object.adata.zarr.zip") })
        ev.register_data_params(adata, { "adata_url": "https://example.com/my_object.adata.zarr" })
        ev.register_data_params(sdata, { "sdata_path": "/path/to/my_object.sdata.zarr" })

    :param obj: AnnData or SpatialData object.
    :type obj: AnnData or SpatialData
    :param dict params: Data-related parameters for the wrapper class.
    """
    if not isinstance(params, dict):
        raise TypeError(f"Expected params to be a dict, got {type(params).__name__} instead.")
    OBJECT_TO_PARAMS[id(obj)] = dict(params)


def _get_param_prefix(obj):
    """
    Get the wrapper parameter prefix for an object.

    :param obj: AnnData or SpatialData object.
    :returns: "adata" for AnnData objects, "sdata" for SpatialData objects.
    """
    from anndata import AnnData
    from spatialdata import SpatialData

    if isinstance(obj, AnnData):
        return "adata"
    if isinstance(obj, SpatialData):
        return "sdata"
    raise TypeError(f"Expected an AnnData or SpatialData object, got {type(obj).__name__} instead.")


def register_data_path(obj, obj_path):
    """
    Register an AnnData or SpatialData object with its corresponding local file path.
    Convenience wrapper around :func:`register_data_params`, equivalent to
    ``register_data_params(adata, { "adata_path": obj_path })``
    (or ``{ "sdata_path": obj_path }`` for SpatialData objects).

    .. code-block:: python

        ev.register_data_path(adata, "/path/to/my_object.adata.zarr")
        ev.register_data_path(adata, "/path/to/my_object.h5ad")
        ev.register_data_path(sdata, "/path/to/my_object.sdata.zarr")

    :param obj: AnnData or SpatialData object.
    :type obj: AnnData or SpatialData
    :param str obj_path: File path where the object is stored.
    """
    register_data_params(obj, {f"{_get_param_prefix(obj)}_path": obj_path})


def register_data_store(obj, obj_store):
    """
    Register an AnnData or SpatialData object with a Zarr store containing the object.
    Convenience wrapper around :func:`register_data_params`, equivalent to
    ``register_data_params(adata, { "adata_store": obj_store })``
    (or ``{ "sdata_store": obj_store }`` for SpatialData objects).

    .. code-block:: python

        ev.register_data_store(adata, zarr.storage.ZipStore("/path/to/my_object.adata.zarr.zip"))

    :param obj: AnnData or SpatialData object.
    :type obj: AnnData or SpatialData
    :param obj_store: Zarr store (or local path to a Zarr store) where the object is stored.
    :type obj_store: zarr.storage.StoreLike or str
    """
    register_data_params(obj, {f"{_get_param_prefix(obj)}_store": obj_store})


def register_data_url(obj, obj_url):
    """
    Register an AnnData or SpatialData object with a remote URL where the object is hosted.
    Convenience wrapper around :func:`register_data_params`, equivalent to
    ``register_data_params(adata, { "adata_url": obj_url })``
    (or ``{ "sdata_url": obj_url }`` for SpatialData objects).
    For AnnData objects stored as H5AD files, the reference spec JSON file is assumed to be
    hosted alongside the H5AD file (i.e., ``my_object.h5ad`` -> ``my_object.ref.json``).

    .. code-block:: python

        ev.register_data_url(adata, "https://example.com/my_object.adata.zarr")
        ev.register_data_url(sdata, "https://example.com/my_object.sdata.zarr")

    :param obj: AnnData or SpatialData object.
    :type obj: AnnData or SpatialData
    :param str obj_url: URL where the object is hosted.
    """
    prefix = _get_param_prefix(obj)
    params = {f"{prefix}_url": obj_url}
    if prefix == "adata" and obj_url.endswith(".h5ad"):
        params["ref_url"] = f"{obj_url[:-len('.h5ad')]}.ref.json"
    register_data_params(obj, params)


def _get_registered_params(obj):
    """
    Get a copy of the data-related parameters registered for an object, if any.

    :param obj: AnnData or SpatialData object.
    :returns: Dictionary of parameters, or None if the object has not been registered.
    """
    if id(obj) in OBJECT_TO_PARAMS:
        return dict(OBJECT_TO_PARAMS[id(obj)])
    return None


def _get_adata_h5ad_ref_json_filepath(h5ad_path):
    """
    Create a Reference Spec JSON filepath for an AnnData object stored in H5AD format.

    :param str h5ad_path: H5AD filepath.
    :returns: Reference Spec JSON filepath.
    """
    base, ext = os.path.splitext(h5ad_path)
    if ext != ".h5ad":
        raise ValueError(f"Expected .h5ad file extension, got {ext} instead.")
    json_filepath = f"{base}.ref.json"

    is_url = config.get('data.wrapper_param_suffix') == '_url'
    if is_url:
        # Assume that the JSON file is hosted alongside the H5AD file correctly.
        return json_filepath

    if not isfile(json_filepath) or config.get('data.overwrite') is True:
        ref_dict = generate_h5ad_ref_spec(h5ad_path)
        with open(json_filepath, "w") as f:
            json.dump(ref_dict, f)
    return json_filepath


def _get_adata_filepath(adata):
    """
    Writes an AnnData object to a new file in the output directory.

    :param AnnData adata: AnnData object.
    :returns: Zarr or H5AD filepath.
    """
    is_url = config.get('data.wrapper_param_suffix') == '_url'
    if is_url:
        raise ValueError("Using AnnData with adata_url requires using register_data_params to register the AnnData object with its corresponding URL, e.g., register_data_params(adata, { 'adata_url': url }).")
        
    # If not, create a new filepath and write the AnnData object to disk.
    # TODO: support Zipped zarr stores?
    is_zarr = config.get('data.anndata_format') == 'zarr'

    file_ext = ".adata.zarr" if is_zarr else ".h5ad"
    file_name = f"{str(uuid4())}{file_ext}"
    out_filepath = join(config.get('data.out_dir'), file_name)

    # Check if the file/folder already exists.
    if exists(out_filepath) and ((is_zarr and isdir(out_filepath)) or (not is_zarr and isfile(out_filepath))):
        if config.get('data.overwrite') is True:
            shutil.rmtree(out_filepath)
        else:
            warnings.warn(f"File {out_filepath} already exists. To overwrite, set data.overwrite to True in the easy_vitessce configuration. By not overwriting, there is a risk of using stale data.")
            # Return early, as the user does not want to overwrite.
            return out_filepath

    os.makedirs(config.get('data.out_dir'), exist_ok=True)
    if is_zarr:
        adata.write_zarr(out_filepath, chunks=[adata.shape[0], VAR_CHUNK_SIZE])
    else:
        adata.write_h5ad(out_filepath)

    return out_filepath

def _get_adata_wrapper_params(adata):
    """
    Get a partial set of parameters for AnnDataWrapper based on the AnnData object and the current configuration.
    
    :param AnnData adata: AnnData object.
    :returns: Dictionary of data-related parameters for AnnDataWrapper.
    """
    registered_params = _get_registered_params(adata)
    if registered_params is not None:
        # If the user registered an H5AD file path without a reference spec,
        # generate the reference spec JSON alongside it.
        adata_path = registered_params.get('adata_path')
        has_ref = 'ref_path' in registered_params or 'ref_url' in registered_params
        if isinstance(adata_path, str) and adata_path.endswith('.h5ad') and not has_ref:
            registered_params['ref_path'] = _get_adata_h5ad_ref_json_filepath(adata_path)
        return registered_params

    # The AnnData object has not been registered, so write it to disk.
    # TODO: do a basic check to determine whether the data on-disk is stale compared to the AnnData object in-memory.
    adata_path = _get_adata_filepath(adata)

    if adata_path.endswith('.zarr'):
        is_store = config.get('data.wrapper_param_suffix') == '_store'
        result = {'adata_store' if is_store else 'adata_path': adata_path}
    else:
        assert adata_path.endswith('.h5ad'), "Expected an H5AD file extension."
        result = {
            'adata_path': adata_path,
            'ref_path': _get_adata_h5ad_ref_json_filepath(adata_path),
        }

    # Register the parameters so that subsequent calls reuse the written file.
    register_data_params(adata, result)
    return dict(result)


# SpatialData analogs of above functions.

def _get_sdata_filepath(sdata):
    """
    Writes a SpatialData object to a new Zarr store in the output directory.

    :param SpatialData sdata: SpatialData object.
    :returns: Zarr filepath.
    """
    # Create a new filepath and write the SpatialData object to disk.
    # TODO: support Zipped zarr stores?
    file_ext = ".sdata.zarr"
    file_name = f"{str(uuid4())}{file_ext}"
    out_filepath = join(config.get('data.out_dir'), file_name)

    # Check if the file/folder already exists.
    if exists(out_filepath) and isdir(out_filepath):
        if config.get('data.overwrite') is True:
            shutil.rmtree(out_filepath)
        else:
            warnings.warn(f"File {out_filepath} already exists. To overwrite, set data.overwrite to True in the easy_vitessce configuration. By not overwriting, there is a risk of using stale data.")
            # Return early, as the user does not want to overwrite.
            return out_filepath

    os.makedirs(config.get('data.out_dir'), exist_ok=True)
    
    sdata.write(out_filepath, overwrite=config.get('data.overwrite'))

    return out_filepath

def _get_sdata_wrapper_params(sdata):
    """
    Get a partial set of parameters for SpatialDataWrapper based on the SpatialData object and the current configuration.

    :param SpatialData sdata: SpatialData object.
    :returns: Dictionary of data-related parameters for SpatialDataWrapper.
    """
    registered_params = _get_registered_params(sdata)
    if registered_params is not None:
        return registered_params

    is_url = config.get('data.wrapper_param_suffix') == '_url'
    if is_url:
        raise ValueError("Using SpatialData with sdata_url requires using register_data_params to register the SpatialData object with its corresponding URL, e.g., register_data_params(sdata, { 'sdata_url': url }).")

    is_store = config.get('data.wrapper_param_suffix') == '_store'
    param_name = 'sdata_store' if is_store else 'sdata_path'

    # If the SpatialData object is already backed and self-contained, use its on-disk location.
    # This is not registered, since sdata.path can change (e.g., after sdata.write).
    if sdata.is_backed() and sdata.is_self_contained():
        return {param_name: str(sdata.path)}

    # Otherwise, write it to disk.
    # TODO: do a basic check to determine whether the data on-disk is stale compared to the SpatialData object in-memory.
    result = {param_name: _get_sdata_filepath(sdata)}

    # Register the parameters so that subsequent calls reuse the written store.
    register_data_params(sdata, result)
    return dict(result)
