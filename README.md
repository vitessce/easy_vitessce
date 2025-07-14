# Easy Vitessce

A package to easily use Vitessce to create interactive plots for single-cell data.

**Supported plots:**

- UMAP
- PCA
- t-SNE
- Spatial
- Violin
- Dotplot
- Heatmap


# Installation

Install package using pip: 

`pip install 'easy_vitessce @ git+https://github.com/luoselena/easy_vitessce@main'`

# How to Use

### Importing Easy Vitessce: 

`from easy_vitessce import configure_plots`

Scanpy is also required to run the package: `import scanpy as sc`

All interactive plots are **enabled by default.**

### Supported Functions and Parameters

<ins> **embedding:** </ins>

`basis:` Plot type. "umap", "pca", or "tsne".

`color:` gene, e.g. "CD79A". 

`color_map:` colormap. "viridis", "plasma", or "jet".

<ins> **spatial**: </ins>

`color:` annotations of observations, e.g. "log1p_n_genes_by_counts". 

`color_map:` colormap. "viridis", "plasma", or "jet".

<ins> **violin**: </ins>


<ins> **dotplot**: </ins>

`markers:` list of genes, e.g.  ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']. 

`groupby:` observation grouping, e.g. "bulk_labels"

<ins> **heatmap**: </ins>

`color_map:` colormap. "viridis", "plasma", or "jet".

`markers:` list of genes.

`groupby:` observation grouping.

### Deactivating Interactive Plots:

`configure_plots(disable_plots = ["spatial", "violin"]`

### Reactivating Interactive Plots:
`configure_plots(enable_plots = ["spatial", "violin"]`

