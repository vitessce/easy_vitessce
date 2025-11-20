# Easy Vitessce

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/vitessce/easy_vitessce/blob/main/docs/notebooks/scanpy_pbmc68k.ipynb)


 🪄 *Configure Vitessce with a single line of code!*
 
Turn your static [Scanpy](https://github.com/scverse/scanpy) and [SpatialData](https://github.com/scverse/spatialdata-plot) plots into interactive [Vitessce](https://github.com/vitessce/vitessce) visualizations simply by importing the `easy_vitessce` package!


**Supported Functions**

- `sc.pl.umap`
- `sc.pl.tsne`
- `sc.pl.pca`
- `sc.pl.diffmap`
- `sc.pl.embedding`
- `sc.pl.violin`
- `sc.pl.dotplot`
- `sc.pl.heatmap`
- `sdata.pl` (`.render_images`, `.render_labels`, `.render_shapes`, `.render_points`)

See the [example notebooks](./docs/notebooks) and the [documentation](https://vitessce.github.io/easy_vitessce/) for further details.


## Installation

Install package using pip: 

```sh
pip install easy_vitessce
```

## How to Use

#### Importing Easy Vitessce

```py
import easy_vitessce as ev
```

🪄 By default, interactive plots are **enabled** via this import statement.

#### Deactivating Interactive Plots:

```py
ev.disable_plots(["embedding", "violin", "spatialdata-plot"])
# or, to disable all interactive plots and return to static plotting mode
ev.disable_plots()
```

#### Reactivating Interactive Plots:

```py
ev.enable_plots(["spatialdata-plot"])
# or, to enable all interactive plots
ev.enable_plots()
```

## Troubleshooting

See the [Troubleshooting](https://github.com/vitessce/vitessce-python?tab=readme-ov-file#troubleshooting) section of the `vitessce-python` repository for tips.

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
open docs/_build/html/index.html
```

### Launch Jupyter notebook or lab

```sh
# uv sync --extra dev
uv run jupyter notebook --notebook-dir .
# or
uv run jupyter lab --notebook-dir .
```

## Citation

To cite EasyVitessce in your work, please use:

```bibtex
@article{luo2025easyvitessce,
  title = {{EasyVitessce: auto-magically adding interactivity to Scverse single-cell and spatial biology plots}},
  author = {Luo, Selena and Keller, Mark S. and Kakar, Tabassum and Choy, Lisa and Gehlenborg, Nils},
  journal = {arXiv},
  year = {2025},
  month = oct,
  doi = {10.48550/arXiv.2510.19532}
}
```

To cite Vitessce in your work, please use:

```bibtex
@article{keller2024vitessce,
  title = {{Vitessce: integrative visualization of multimodal and spatially resolved single-cell data}},
  author = {Keller, Mark S. and Gold, Ilan and McCallum, Chuck and Manz, Trevor and Kharchenko, Peter V. and Gehlenborg, Nils},
  journal = {Nature Methods},
  year = {2024},
  month = sep,
  doi = {10.1038/s41592-024-02436-x}
}
```

If you use the image rendering functionality, please additionally cite Viv:

```bibtex
@article{manz2022viv,
  title = {{Viv: multiscale visualization of high-resolution multiplexed bioimaging data on the web}},
  author = {Manz, Trevor and Gold, Ilan and Patterson, Nathan Heath and McCallum, Chuck and Keller, Mark S. and Herr, II, Bruce W. and Börner, Kay and Spraggins, Jeffrey M. and Gehlenborg, Nils},
  journal = {Nature Methods},
  year = {2022},
  month = may,
  doi = {10.1038/s41592-022-01482-7}
}
```
