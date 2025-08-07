# Easy Vitessce

 🪄 *Magic in a single line of code!* 
 
Turn your static [Scanpy](https://github.com/scverse/scanpy) and [SpatialData](https://github.com/scverse/spatialdata) plots into interactive [Vitessce](https://github.com/vitessce/vitessce) visualizations with _Easy Vitessce_ for spatial and single-cell data just by adding the `easy_vitessce` package!

**Supported Plots**

- UMAP
- PCA
- t-SNE
- Spatial (Scanpy version)
- Spatial (SpatialData version)
- Violin
- Dotplot
- Heatmap


## Installation

Install package using pip: 

```
pip install 'easy_vitessce @ git+https://github.com/vitessce/easy_vitessce@main'
```

## How to Use

#### Importing Easy Vitessce

```
from easy_vitessce import configure_plots
```

🪄 All interactive plots are **enabled magically**.

_Note that Scanpy is also required to run the package:_ 
```
import scanpy as sc
```

#### Supported Functions and Parameters

<ins> **embedding:** </ins>

`basis:` Plot type. "umap", "pca", or "tsne".

`color:` gene, e.g. "CD79A". 

`color_map:` color map. "viridis", "plasma", or "jet".

<ins> **spatial (Scanpy ver.)**: </ins>

`color:` annotations of observations, e.g. "log1p_n_genes_by_counts". 

`color_map:` color map. "viridis", "plasma", or "jet".

<ins> **Spatial (SpatialData ver.)** </ins>

`spatialdata_filepath:` filepath of spatialdata zarr file containing image data.

`zip_filepath:` filepath of zip folder.

`render_images():` renders image.
* `element:` image data location.

`render_shapes():` renders shapes, e.g. spots.
* `element:` element to be rendered, e.g. "cells".
* `color:` gene.
* `color_map:` color map. "viridis", "plasma", or "jet".

`render_labels():` renders labels.
* `element:` label data location.

`show():` displays interactive plot.

<ins> **violin**: </ins>

`markers:` list of genes.
`groupby:` observation grouping, e.g. "bulk_labels"

<ins> **dotplot**: </ins>

`markers:` list of genes, e.g.  ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']. 

`groupby:` observation grouping, e.g. "bulk_labels"

<ins> **heatmap**: </ins>

`color_map:` colormap. "viridis", "plasma", or "jet".

`markers:` list of genes.

`groupby:` observation grouping.

#### Deactivating Interactive Plots:

`configure_plots(disable_plots = ["spatial", "violin"]`

#### Reactivating Interactive Plots:
`configure_plots(enable_plots = ["spatial", "violin"]`

