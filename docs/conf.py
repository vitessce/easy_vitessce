# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

#sys.path.insert(0, os.path.abspath('..'))
#sys.path.append("/easy_vitessce/src/easy_vitessce")

import easy_vitessce
from easy_vitessce import configure_plots
from easy_vitessce.spatialdata_plot import VitesscePlotAccessor

project = 'Easy Vitessce'
copyright = '2025, HIDIVE Lab'
author = 'HIDIVE Lab'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.todo", "sphinx.ext.viewcode", "sphinx.ext.autodoc", 'sphinx_copybutton', 'sphinx.ext.githubpages', 'myst_parser', 'nbsphinx']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'notebooks/data/**']

source_suffix = {
  '.rst': 'restructuredtext',
  '.md': 'markdown',
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']
html_theme_options = {
  "show_toc_level": 3,
  "max_navbar_depth": 2
}

# -- Options for nbsphinx -------------------------------------------------

# By default, notebooks are rendered as-is (without executing them).
nbsphinx_execute = 'never'
#nbsphinx_allow_errors = True

# Notebooks ending in `_remote.ipynb` reference remotely-hosted data
# (via `ev.register_data_params`), so their interactive widget outputs
# can function in the static HTML site. These notebooks are executed
# during the docs build (e.g., in CI) so that their outputs
# (including the Vitessce widget state) are embedded in the built HTML.
NBSPHINX_EXECUTE_SUFFIX = '_remote.ipynb'

# Allow time for downloading example datasets.
nbsphinx_timeout = 1800


def _execute_remote_notebooks(app, docname, source):
    import json

    srcfile = str(app.env.doc2path(docname))
    if not srcfile.endswith(NBSPHINX_EXECUTE_SUFFIX):
        return
    nb = json.loads(source[0])
    # Per-notebook override of `nbsphinx_execute`.
    # Reference: https://nbsphinx.readthedocs.io/en/latest/never-execute.html
    nb.setdefault('metadata', {}).setdefault('nbsphinx', {})['execute'] = 'always'
    source[0] = json.dumps(nb)


def setup(app):
    # The source-read event provides the raw notebook JSON before nbsphinx parses it.
    app.connect('source-read', _execute_remote_notebooks)