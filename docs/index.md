# Easy Vitessce Documentation

Easy Vitessce is a Python package for turning Scanpy and SpatialData plots into interactive [Vitessce](https://vitessce.io/) visualizations with minimal code changes.


## Installation
Requires Python 3.10 or greater.

```sh
pip install easy_vitessce
```


## How to Use and Examples

The package can be imported with

```py
import easy_vitessce as ev
```

By default, interactive plots are enabled via running this import statement.

### Deactivating/Reactivating Interactive Plots

Run `ev.disable_plots` to deactivate Vitessce plots.

Run `ev.enable_plots` to re-activate Vitessce plots.

```py
ev.disable_plots(["embedding", "violin", "spatialdata-plot"])

ev.enable_plots(["spatialdata-plot", "violin"])
```


### Spatial Plot (SpatialData-Plot version)

**Note:** This example uses SpatialData's [mouse brain MERFISH dataset.](https://spatialdata.scverse.org/en/latest/tutorials/notebooks/datasets/README.html)

```py
sdata = sd.read_zarr(spatialdata_filepath)
sdata.pl.render_images(element="rasterized").pl.render_shapes(element="cells", color="Acta2").pl.show()
```

`spatialdata_filepath` should lead to a `.zarr` file containing spatial data with an `Images` folder. The file structure of the example above is as follows. Since it does not have a `Labels` folder, calling `pl.render_labels()` will not display any data.

```
SpatialData object, with associated Zarr store:
├── Images
│     └── 'rasterized': DataArray[cyx] (1, 522, 575)
├── Points
│     └── 'single_molecule': DataFrame with shape: (<Delayed>, 3) (2D points)
├── Shapes
│     ├── 'anatomical': GeoDataFrame shape: (6, 1) (2D shapes)
│     └── 'cells': GeoDataFrame shape: (2389, 2) (2D shapes)
└── Tables
      └── 'table': AnnData (2389, 268)
```

<p>
    <img alt="static_sd" src="_static/static_sd_spatial.png" style="width: 30%" />
    <img alt="right_arrow1" src="_static/right_arrow_transparent.png" style="width: 8%" />
    <img alt="vitessce_sd" src="_static/spatial_documentation.gif" style="width: 60%" />
</p>

### Spatial Plot (Scanpy version)

Easy Vitessce's `spatial` function also displays a spatial plot, but with Scanpy's syntax. This example uses Scanpy's [Visium dataset.](https://scanpy.readthedocs.io/en/stable/generated/scanpy.datasets.visium_sge.html)

```py
adata = sc.datasets.visium_sge(sample_id="Targeted_Visium_Human_Glioblastoma_Pan_Cancer", include_hires_tiff=True)

sc.pl.spatial(adata, color = "log1p_n_genes_by_counts")
```

<p>
    <img alt="static_sc_spatial" src="_static/static_sc_spatial.png" style="width: 25%" />
    <img alt="right_arrow2" src="_static/right_arrow_transparent.png" style="width: 10%" />
    <img alt="sc_spatial_documentation" src="_static/sc_spatial_documentation.gif" style="width: 60%" />
</p>

### Scatterplots

Easy Vitessce's `embedding` function displays UMAP, PCA, and t-SNE scatterplots.

The functions `umap`, `pca`, and `tsne` can also be used for convenience.

```py
adata = sc.datasets.pbmc68k_reduced()

sc.pl.embedding(adata, basis="X_umap", color="CD79A")
sc.pl.embedding(adata, basis="X_pca", color=["CD79A", "CD53"])
sc.pl.embedding(adata, basis="X_tsne", color=["bulk_labels", "louvain", "phase"])

sc.pl.umap(adata, color="CD79A")
sc.pl.pca(adata, color="CD79A")
sc.pl.tsne(adata, color="CD79A")
```

<p>
    <img alt="static_umap" src="_static/static_umap.png" style="width: 37%" />
    <img alt="right_arrow3" src="_static/right_arrow_transparent.png" style="width: 6%" />
    <img alt="vitessce_umap" src="_static/updated_umap.gif" style="width: 55%" />
</p>


### Dotplot

**Note:** To select/deselect multiple genes, hold SHIFT while clicking on genes in the Gene List.

```py
adata = sc.datasets.pbmc68k_reduced()

sc.pl.dotplot(adata, var_names=["C1QA", "PSAP", "CD79A", "CD79B", "CST3", "LYZ"], groupby="bulk_labels")
```

<p>
    <img alt="static_dotplot" src="_static/static_dotplot.png" style="width: 37%" />
    <img alt="right_arrow4" src="_static/right_arrow_transparent.png" style="width: 6%" />
    <img alt="vitessce_dotplot" src="_static/dotplot_example.gif" style="width: 55%" />
</p>

### Violin Plot

```py
adata = sc.datasets.pbmc68k_reduced()

sc.pl.violin(adata, keys="AP2S1", groupby="bulk_labels")
```

<p>
    <img alt="static_violin" src="_static/static_violin.png" style="width: 37%" />
    <img alt="right_arrow5" src="_static/right_arrow_transparent.png" style="width: 6%" />
    <img alt="vitessce_violin" src="_static/updated_violin.gif" style="width: 55%" />
</p>

### Heatmap

```py
adata = sc.datasets.pbmc68k_reduced()

sc.pl.heatmap(adata, groupby="bulk_labels", var_names=['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ'])
```

<p>
    <img alt="static_heatmap" src="_static/static_heatmap.png" style="width: 37%" />
    <img alt="right_arrow6" src="_static/right_arrow_transparent.png" style="width: 6%" />
    <img alt="vitessce_heatmap" src="_static/heatmap.gif" style="width: 55%" />
</p>

<!-- Reference: https://myst-parser.readthedocs.io/en/latest/develop/background.html#the-relationship-between-myst-restructuredtext-and-sphinx -->

```{toctree}
:maxdepth: 2
:hidden:
:caption: Contents

easy_vitessce
examples
customization
advanced
Example notebooks <https://github.com/vitessce/easy_vitessce/tree/main/docs/notebooks>
View on GitHub <https://github.com/vitessce/easy_vitessce>
```

## Citation

To cite EasyVitessce in your work, please use:

> S Luo, MS Keller, T Kakar, L Choy, N Gehlenborg. "EasyVitessce: auto-magically adding interactivity to Scverse single-cell and spatial biology plots", _arXiv_ (2025). doi:[10.48550/arXiv.2510.19532](    
https://doi.org/10.48550/arXiv.2510.19532)