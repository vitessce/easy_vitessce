import os
from os.path import join, isfile, isdir, exists
from uuid import uuid4
import shutil
import warnings
import json

from vitessce.data_utils import VAR_CHUNK_SIZE, generate_h5ad_ref_spec

from .widget import config

OBJECT_TO_PATH = {}

def register_data_path(obj, obj_path):
    """
    Register an AnnData or SpatialData object with its corresponding file path.
    This allows to avoid writing a large AnnData to disk if it has already been written.
    
    :param obj: AnnData object.
    :type obj: AnnData or SpatialData
    :param str obj_path: File path where the object is stored.
    """
    # AnnData objects are not hashable,
    # so we use the id() function to get a unique identifier for the object.
    # Reference: https://github.com/scverse/anndata/issues/742
    OBJECT_TO_PATH[id(obj)] = obj_path


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
    Creates Zarr filepath for AnnData object. Prevents creating multiple files with the same name.

    :param AnnData adata: AnnData object.
    :returns: Zarr filepath.
    """
    # Check if the AnnData object has already been written to disk, by checking if it exists in the dictionary.
    if id(adata) in OBJECT_TO_PATH:
        # TODO: do a basic check to determine whether the data on-disk is stale compared to the AnnData object in-memory.
        return OBJECT_TO_PATH[id(adata)]
    
    is_url = config.get('data.wrapper_param_suffix') == '_url'
    if is_url:
        raise ValueError("Using AnnData with adata_url requires using register_data_path to register the AnnData object with its corresponding URL.")
        
    # If not, create a new filepath and write the AnnData object to disk.
    # TODO: support Zipped zarr stores?
    is_zarr = config.get('data.anndata_format') == 'zarr'

    file_ext = ".adata.zarr" if is_zarr else ".h5ad"
    file_name = f"{str(uuid4())}{file_ext}"
    out_filepath = join(config.get('data.out_dir'), file_name)

    # Add the new AnnData object and its filepath to the weakref dictionary.
    register_data_path(adata, out_filepath)

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
    result = {}

    adata_path = _get_adata_filepath(adata)
    is_zarr = adata_path.endswith('.zarr')
    is_url = adata_path.startswith('http://') or adata_path.startswith('https://')

    if is_zarr:
        is_store = config.get('data.wrapper_param_suffix') == '_store'
        if is_url:
            assert adata_path.startswith(('http://', 'https://')), "Expected a valid URL."
            result['adata_url'] = adata_path
        elif is_store:
            result['adata_store'] = adata_path
        else:
            result['adata_path'] = adata_path
    else:
        assert adata_path.endswith('.h5ad'), "Expected an H5AD file extension."
        json_path = _get_adata_h5ad_ref_json_filepath(adata_path)

        if is_url:
            assert adata_path.startswith(('http://', 'https://')), "Expected a valid URL."
            assert json_path.startswith(('http://', 'https://')), "Expected a valid URL."
            result['adata_url'] = adata_path
            result['ref_url'] = json_path
        else:
            result['adata_path'] = adata_path
            result['ref_path'] = json_path
    
    return result


# SpatialData analogs of above functions.

def _get_sdata_filepath(sdata):
    """
    Creates Zarr filepath for AnnData object. Prevents creating multiple files with the same name.

    :param AnnData adata: AnnData object.
    :returns: Zarr filepath.
    """
    # Check if the AnnData object has already been written to disk, by checking if it exists in the dictionary.
    if id(sdata) in OBJECT_TO_PATH:
        # TODO: do a basic check to determine whether the data on-disk is stale compared to the AnnData object in-memory.
        return OBJECT_TO_PATH[id(sdata)]
    
    is_url = config.get('data.wrapper_param_suffix') == '_url'
    if is_url:
        raise ValueError("Using SpatialData with sdata_url requires using register_data_path to register the SpatialData object with its corresponding URL.")

    # Check if the SpatialData object is already backed and self-contained.
    if sdata.is_backed() and sdata.is_self_contained():
        return str(sdata.path)

    # If not, create a new filepath and write the SpatialData object to disk.
    # TODO: support Zipped zarr stores?
    file_ext = ".sdata.zarr"
    file_name = f"{str(uuid4())}{file_ext}"
    out_filepath = join(config.get('data.out_dir'), file_name)

    # Add the new SpatialData object and its filepath to the weakref dictionary.
    register_data_path(sdata, out_filepath)

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
    result = {}

    sdata_path = _get_sdata_filepath(sdata)
    is_url = sdata_path.startswith('http://') or sdata_path.startswith('https://')

    is_store = config.get('data.wrapper_param_suffix') == '_store'
    if is_url:
        assert sdata_path.startswith(('http://', 'https://')), "Expected a valid URL."
        result['sdata_url'] = sdata_path
    elif is_store:
        result['sdata_store'] = sdata_path
    else:
        result['sdata_path'] = sdata_path

    return result