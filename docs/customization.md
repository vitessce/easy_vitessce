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

To avoid redundantly writing AnnData objects, you can tell `easy_vitessce` that an AnnData object has already been stored on-disk:

```py
import easy_vitessce as ev

ev.register_data_path(adata, '/path/to/my_object.adata.zarr')
# or
ev.register_data_path(adata, '/path/to/my_object.h5ad')
```

For SpatialData objects, this problem is resolved by the fact that disk-backed objects have an `sdata.path` attribute to specify the on-disk location.

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
