# Customization

While Easy Vitessce is designed to work with minimal code changes, and provide sensible defaults, there are some additional ways to customize the behavior.

The simplest form of configuration is to enable or disable plots using `configure_plots`:

```python
import easy_vitessce as ev

# Enable or disable particular plotting functions
ev.configure_plots(enable_plots=["embedding"], disable_plots=["violin", "heatmap"])
# or
ev.enable_plots(["embedding"])
# or
ev.disable_plots(["violin", "heatmap"])
```

For all other types of configuration, import the [Donfig-based](https://donfig.readthedocs.io/en/latest/index.html) `config` variable:

```py
import easy_vitessce as ev

# Pretty-print the current configuration
ev.config.pprint()
```

## Widget Configuration

Configure the parameters that `easy_vitessce` plotting functions will internally pass to the [`VitessceConfig.widget`](https://python-docs.vitessce.io/api_config.html#vitessce-widget) function:

```py
import easy_vitessce as ev

ev.config.update({ 'widget': { 'js_dev_mode': True } })
# or
ev.config.set({ 'widget': { 'js_dev_mode': True } })
```

Configure whether to use `VitessceConfig.widget` or `VitessceConfig.display` to show the widget:

```py
import easy_vitessce as ev

ev.config.update({ 'widget_function': 'display' })
```

## Data-Related Configuration

By default, Easy Vitessce writes AnnData objects to disk using `AnnData.write_zarr` (as there is no way to know whether/where an AnnData object has been stored on-disk, given only an `adata` variable).
If the AnnData object is large, this can be time- and disk-consuming.

To avoid redundantly writing AnnData objects, you can tell `easy_vitessce` where an AnnData object has already been stored, by registering the data-related parameters that should be passed to [AnnDataWrapper](https://python-docs.vitessce.io/api_data.html#vitessce.wrappers.AnnDataWrapper):

```py
import easy_vitessce as ev

ev.register_data_params(adata, { "adata_path": "/path/to/my_object.adata.zarr" })
# or, for H5AD files (a reference spec JSON file will be generated alongside the H5AD file if ref_path is not provided)
ev.register_data_params(adata, { "adata_path": "/path/to/my_object.h5ad" })
```

Because an entire dictionary of parameters is registered, alternative storage mechanisms can be used, such as a Zarr store or a remote URL:

```py
import zarr
import easy_vitessce as ev

ev.register_data_params(adata, { "adata_store": zarr.storage.ZipStore("/path/to/my_object.adata.zarr.zip") })
# or
ev.register_data_params(adata, { "adata_url": "https://example.com/my_object.adata.zarr" })
```

For SpatialData objects, disk-backed objects have an `sdata.path` attribute to specify the on-disk location.
SpatialData objects can also be registered, with parameters for [SpatialDataWrapper](https://python-docs.vitessce.io/api_data.html#vitessce.wrappers.SpatialDataWrapper):

```py
ev.register_data_params(sdata, { "sdata_url": "https://example.com/my_object.sdata.zarr" })
```

The following convenience functions are equivalent to calling `ev.register_data_params` with a single `*_path`, `*_store`, or `*_url` parameter (using the `adata_` or `sdata_` prefix based on the object type):

```py
import easy_vitessce as ev

ev.register_data_path(adata, "/path/to/my_object.adata.zarr")
# Equivalent to ev.register_data_params(adata, { "adata_path": "/path/to/my_object.adata.zarr" })

ev.register_data_store(adata, zarr.storage.ZipStore("/path/to/my_object.adata.zarr.zip"))
# Equivalent to ev.register_data_params(adata, { "adata_store": zarr.storage.ZipStore("/path/to/my_object.adata.zarr.zip") })

ev.register_data_url(sdata, "https://example.com/my_object.sdata.zarr")
# Equivalent to ev.register_data_params(sdata, { "sdata_url": "https://example.com/my_object.sdata.zarr" })
```

For AnnData objects hosted as H5AD files, `ev.register_data_url` assumes the reference spec JSON file is hosted alongside the H5AD file (e.g., `my_object.h5ad` and `my_object.ref.json`).

Other data-related options include configuring usage of `AnnData.write_zarr` (the default) vs. `AnnData.write_h5ad`:

```py
import easy_vitessce as ev

ev.config.update({ 'data.anndata_format': 'h5ad' })
```

Or, specify an alternative output directory (the default is `data`):

```py
import easy_vitessce as ev

ev.config.update({ 'data.out_dir': '/path/to/some/dir' })
```
