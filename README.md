# Easy Vitessce

 🪄 *Configure Vitessce with a single line of code!*
 
Turn your static [Scanpy](https://github.com/scverse/scanpy) and [SpatialData](https://github.com/scverse/spatialdata) plots into interactive [Vitessce](https://github.com/vitessce/vitessce) visualizations with _Easy Vitessce_ for spatial and single-cell data just by adding the `easy_vitessce` package!

🚧 work in progress 🚧

**Supported Functions**

- `sc.pl.umap`
- `sc.pl.tsne`
- `sc.pl.pca`
- `sc.pl.embedding`
- `sc.pl.spatial`
- `sc.pl.violin`
- `sc.pl.dotplot`
- `sc.pl.heatmap`
- `sdata.pl.render_images`
- `sdata.pl.render_labels`
- `sdata.pl.render_shapes`

See the [documentation](https://vitessce.github.io/easy_vitessce/) for further details.


## Installation

Install package using pip: 

```sh
pip install easy_vitessce
```

## How to Use

#### Importing Easy Vitessce

```py
from easy_vitessce import configure_plots
```

🪄 By default, interactive plots are **enabled**.

#### Deactivating Interactive Plots:

```py
configure_plots(disable_plots = ["spatial", "violin"])
```

#### Reactivating Interactive Plots:

```py
configure_plots(enable_plots = ["spatial", "violin"])
```

## Development

### Set up environment

```sh
uv sync --extra dev --extra docs
```

This command should also be run after updating dependencies in `pyproject.toml`.

### Run tests

```sh
# uv sync --extra dev
uv run pytest
```

### Make documentation

```sh
uv run make html # on mac/linux
# uv run make.bat html # on windows
open _build/html/index.html
```

### Launch Jupyter notebook or lab

```sh
# uv sync --extra dev
uv run jupyter notebook --notebook-dir .
# or
uv run jupyter lab --notebook-dir .
```