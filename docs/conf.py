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

nbsphinx_execute = 'never'
#nbsphinx_allow_errors = True