.. Easy Vitessce Documentation documentation master file, created by
   sphinx-quickstart on Fri Jul 25 15:27:51 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Easy Vitessce Documentation
=========================================

Easy Vitessce is a Python package for turning Scanpy and SpatialData plots into interactive `Vitessce <https://vitessce.io/>`_ visualizations with minimal code changes.

Installation and Use
########################
Requires Python 3.9 or greater.

.. code-block:: 

    pip install 'easy_vitessce @ git+https://github.com/luoselena/easy_vitessce@main'

How to Use and Examples
############################

The package can be imported with

.. code-block:: 

   from easy_vitessce.VitessceSpatialData import VitessceSpatialData
   from easy_vitessce import configure_plots

**Note:** All example datasets are from Scanpy. 

Deactivating/Reactivating Interactive Plots
**************************************************
Passing ``disable_plots`` into ``configure_plots`` will deactivate Vitessce plots.

Passing ``enable_plots`` into ``configure_plots`` will reactivate Vitessce plots.

**Note:** While ``diable_plots`` and ``enable_plots`` can be passed in at the same time, listing the same plot in both will result in an error.

.. code-block:: 

   configure_plots(disable_plots = ["spatial", "violin"])
   
   configure_plots(enable_plots = ["spatial", "violin"])
  
Spatial Plot (SpatialData version)
***********************************************

.. code-block:: 

   spatial_plot = VitessceSpatialData(spatialdata_filepath)
   spatial_plot.pl.render_images(element="rasterized").pl.render_shapes(element="cells", color="Acta2").pl.show()

``spatialdata_filepath`` should contain spatial data with an ``Images`` folder. The file structure of the example above is as follows. Since it does not have a ``Labels`` folder, calling ``pl.render_labels()`` will not display any data.

.. code-block:: 

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

Spatial Plot (Scanpy version)
********************************************
Easy Vitessce's ``spatial`` function also displays a spatial plot, but with Scanpy's syntax.

.. code-block:: 

   adata = sc.datasets.visium_sge(sample_id="V1_Human_Heart", include_hires_tiff=True)

   sc.pl.spatial(adata, color = "log1p_n_genes_by_counts")

Scatterplots
**************************
Easy Vitessce's ``embedding`` function displays UMAP, PCA, and t-SNE scatterplots.

.. code-block:: 

   adata = sc.datasets.pbmc68k_reduced()

   sc.pl.embedding(adata, basis="umap", color="CD79A")
   sc.pl.embedding(adata, basis="pca", color=["CD79A", "CD53"])
   sc.pl.embedding(adata, basis="tsne", color=["bulk_labels", "louvain", "phase"])

Dotplot
***********************

.. code-block:: 

   adata = sc.datasets.pbmc68k_reduced()

   sc.pl.dotplot(adata, markers = ["C1QA", "PSAP", "CD79A", "CD79B", "CST3", "LYZ"], groupby="bulk_labels")

Violin Plot
**************************

.. code-block:: 

   adata = sc.datasets.pbmc68k_reduced()

   sc.pl.violin(adata, markers = "ABI1")

Heatmap
*********************

.. code-block:: 

   adata = sc.datasets.pbmc68k_reduced()

   sc.pl.heatmap(adata, groupby = "bulk_labels", markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ'])


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules

