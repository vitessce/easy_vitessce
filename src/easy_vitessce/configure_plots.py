

from os.path import join
import warnings

import pandas as pd
import numpy as np
import scanpy as sc

from vitessce import (
    VitessceConfig,
    AnnDataWrapper,
    ImageOmeZarrWrapper,
    CoordinationLevel as CL,
    ViewType as vt,
    vconcat,
    hconcat
)

from vitessce.data_utils import (
    VAR_CHUNK_SIZE,
    rgb_img_to_ome_zarr
)

from anndata import AnnData
import spatialdata as sd
from spatialdata import SpatialData
from xarray.core.extensions import _CachedAccessor

from easy_vitessce.spatialdata_plot import VitesscePlotAccessor
from easy_vitessce.widget import _to_widget
from easy_vitessce.data import _get_adata_wrapper_params


def umap(adata, **kwargs):
  """
  Creates interactive UMAP plot.

  :param AnnData adata: AnnData object.
  :param str basis: Name of embedding basis to use (e.g., umap, pca, or tsne). Will look up coordinates in adata.obsm["X_{basis}"].
  :param str color: Gene ID from adata.var.index or categorical label from adata.obs.columns.
  :param str gene_symbols: Key in adata.var to use for gene symbols, if different from adata.var.index.
  :param str layer: Layer in AnnData to use for expression values. Defaults to None, which uses adata.X.
  :param str color_map: Color map (viridis, plasma, jet). Defaults to viridis.
  :param size: Size of dots.
  :type size: float or int
  :param float vmin: Minimum value for color map scaling.
  :param float vmax: Maximum value for color map scaling.
  :param title: Title of the plot.
  :type title: str or list[str]
  :param int ncols: Number of columns to use when laying out multiple scatterplots (when multiple colors are specified). Defaults to 4.
  :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
  """
  return embedding(adata, basis="umap", **kwargs)

def tsne(adata, **kwargs):
  """
  Creates interactive t-SNE plot.

  :param AnnData adata: AnnData object.
  :param str basis: Name of embedding basis to use (e.g., umap, pca, or tsne). Will look up coordinates in adata.obsm["X_{basis}"].
  :param str color: Gene ID from adata.var.index or categorical label from adata.obs.columns.
  :param str gene_symbols: Key in adata.var to use for gene symbols, if different from adata.var.index.
  :param str layer: Layer in AnnData to use for expression values. Defaults to None, which uses adata.X.
  :param str color_map: Color map (viridis, plasma, jet). Defaults to viridis.
  :param size: Size of dots.
  :type size: float or int
  :param float vmin: Minimum value for color map scaling.
  :param float vmax: Maximum value for color map scaling.
  :param title: Title of the plot.
  :type title: str or list[str]
  :param int ncols: Number of columns to use when laying out multiple scatterplots (when multiple colors are specified). Defaults to 4.
  :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
  """
  return embedding(adata, basis="tsne", **kwargs)

def pca(adata, **kwargs):
  """
  Creates interactive PCA plot.

  :param AnnData adata: AnnData object.
  :param str basis: Name of embedding basis to use (e.g., umap, pca, or tsne). Will look up coordinates in adata.obsm["X_{basis}"].
  :param str color: Gene ID from adata.var.index or categorical label from adata.obs.columns.
  :param str gene_symbols: Key in adata.var to use for gene symbols, if different from adata.var.index.
  :param str layer: Layer in AnnData to use for expression values. Defaults to None, which uses adata.X.
  :param str color_map: Color map (viridis, plasma, jet). Defaults to viridis.
  :param size: Size of dots.
  :type size: float or int
  :param float vmin: Minimum value for color map scaling.
  :param float vmax: Maximum value for color map scaling.
  :param title: Title of the plot.
  :type title: str or list[str]
  :param int ncols: Number of columns to use when laying out multiple scatterplots (when multiple colors are specified). Defaults to 4.
  :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
  """
  return embedding(adata, basis="pca", **kwargs)

def diffmap(adata, **kwargs):
  """
  Creates interactive Diffusion Map plot.

  :param AnnData adata: AnnData object.
  :param str basis: Name of embedding basis to use (e.g., umap, pca, or tsne). Will look up coordinates in adata.obsm["X_{basis}"].
  :param str color: Gene ID from adata.var.index or categorical label from adata.obs.columns.
  :param str gene_symbols: Key in adata.var to use for gene symbols, if different from adata.var.index.
  :param str layer: Layer in AnnData to use for expression values. Defaults to None, which uses adata.X.
  :param str color_map: Color map (viridis, plasma, jet). Defaults to viridis.
  :param size: Size of dots.
  :type size: float or int
  :param float vmin: Minimum value for color map scaling.
  :param float vmax: Maximum value for color map scaling.
  :param title: Title of the plot.
  :type title: str or list[str]
  :param int ncols: Number of columns to use when laying out multiple scatterplots (when multiple colors are specified). Defaults to 4.
  :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
  """
  return embedding(adata, basis="diffmap", **kwargs)

# Reference:
# - https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.embedding.html
# - https://github.com/scverse/scanpy/blob/cf8b46dea735c35a629abfaa2e1bab9047281e34/src/scanpy/plotting/_tools/scatterplots.py#L70
def embedding(
        adata,
        basis,
        *,
        color = None,
        gene_symbols = None,
        layer = None,
        color_map = None,
        size = None,
        vmin = None,
        vmax = None,
        title = None,
        ncols = 4,
        **kwargs
    ):
    """
    Creates interactive versions of UMAP, PCA, t-SNE plots.

    :param AnnData adata: AnnData object.
    :param str basis: Name of embedding basis to use (e.g., umap, pca, or tsne). Will look up coordinates in adata.obsm["X_{basis}"].
    :param str color: Gene ID from adata.var.index or categorical label from adata.obs.columns.
    :param str gene_symbols: Key in adata.var to use for gene symbols, if different from adata.var.index.
    :param str layer: Layer in AnnData to use for expression values. Defaults to None, which uses adata.X.
    :param str color_map: Color map (viridis, plasma, jet). Defaults to viridis.
    :param size: Size of dots.
    :type size: float or int
    :param float vmin: Minimum value for color map scaling.
    :param float vmax: Maximum value for color map scaling.
    :param title: Title of the plot.
    :type title: str or list[str]
    :param int ncols: Number of columns to use when laying out multiple scatterplots (when multiple colors are specified). Defaults to 4.
    :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
    """
    basis_name = "t-SNE" if basis == "tsne" else basis.upper()
    basis_obsm_key = f"X_{basis}"

    if basis_obsm_key not in adata.obsm:
        raise ValueError(f"{basis_name} coordinates not found in adata.obsm.")

    # Normalize the color parameter to be either list or None.
    if type(color) == str:
        color = [color]

    if size is None:
        # TODO: make this dependent on number of cells?
        size = 2.5
    
    vc = VitessceConfig(schema_version="1.0.18", name=f"sc.pl.embedding for {basis_name}")

    wrapper_params = {
        'obs_embedding_paths': [f"obsm/{basis_obsm_key}"],
        'obs_embedding_names': [basis_name],
    }

    # Configure wrapper_params.
    has_color_in_var = False
    has_color_in_obs = False

    obs_set_paths = []
    obs_set_names = []
    
    if color is not None:
        for c in color:
            if c in adata.obs.columns:
                has_color_in_obs = True
                obs_set_paths.append(f"obs/{c}")
                obs_set_names.append(c.capitalize())
            elif c in adata.var.index:
                has_color_in_var = True

        if len(obs_set_paths) > 0:
            wrapper_params['obs_set_paths'] = obs_set_paths
            wrapper_params['obs_set_names'] = obs_set_names
        
        if has_color_in_var:
            wrapper_params['obs_feature_matrix_path'] = "X" if layer is None else f"layers/{layer}"
    
    if gene_symbols is not None:
        wrapper_params['feature_labels_path'] = f"var/{gene_symbols}"
    
    dataset = vc.add_dataset(name='embedding data').add_object(AnnDataWrapper(
        **_get_adata_wrapper_params(adata),
        **wrapper_params
    ))

    if type(color) == list and len(color) > 1:
        # Multiple colors: multiple scatterplot views.

        scatterplot_views = []
        obs_sets_view = None
        feature_list_view = None

        for color_i, color_val in enumerate(color):
            scatterplot_props = {}
            if type(title) == str and color_i == 0:
                scatterplot_props['title'] = title
            elif type(title) == list and len(color) == len(title):
                scatterplot_props['title'] = title[color_i]
            
            scatterplot_view = vc.add_view(vt.SCATTERPLOT, dataset=dataset).set_props(**scatterplot_props)
            scatterplot_views.append(scatterplot_view)

            if type(size) is list and len(size) == len(color):
                vc.link_views(scatterplot_views,
                    ["embeddingObsRadiusMode", "embeddingObsRadius"],
                    ["manual", size[color_i]]
                )
            
            # Check if color[i] is in obs or var.
            if color_val in adata.obs.columns:
                if obs_sets_view is None:
                    obs_sets_view = vc.add_view(vt.OBS_SETS, dataset=dataset)
                
                # TODO: simplify this once https://github.com/vitessce/vitessce/issues/2254 is addressed.
                group_name = color_val.capitalize()
                obs_set_categories = adata.obs[color_val].unique().tolist()
                obs_set_paths = [
                    [group_name, str(cat)] for cat in obs_set_categories
                ]

                vc.link_views([scatterplot_view, obs_sets_view],
                    ["obsColorEncoding", "obsSetSelection", "obsSetExpansion"],
                    ["cellSetSelection", obs_set_paths, [[group_name]]]
                )
            elif color_val in adata.var.index:
                if feature_list_view is None:
                    feature_list_view = vc.add_view(vt.FEATURE_LIST, dataset=dataset)
                vc.link_views([scatterplot_view, feature_list_view],
                    ["featureSelection", "obsColorEncoding"],
                    [[color_val], "geneSelection"]
                )
        
        # Link the zoom/pan interactions across all scatterplot views.
        vc.link_views(scatterplot_views,
            ["embeddingZoom", "embeddingTargetX", "embeddingTargetY"],
            [None, None, None]
        )

        # Other linkings that are common across all scatterplot views.
        vc.link_views(scatterplot_views,
            ["embeddingType"],
            [basis_name]
        )

        if size is not None and type(size) is not list:
            vc.link_views(scatterplot_views,
                ["embeddingObsRadiusMode", "embeddingObsRadius"],
                ["manual", size]
            )
        
        if vmin is not None or vmax is not None:
            vc.link_views(scatterplot_views,
                ["featureValueColormapRange"],
                [
                    [
                        vmin if vmin is not None else 0.0,
                        vmax if vmax is not None else 1.0,
                    ]
                ]
            )

        if color_map is not None and color_map in ["viridis", "plasma", "jet", "greys"]:
            vc.link_views(scatterplot_views,
                ["featureValueColormap"],
                [color_map]
            )
        
        # Layout the views.
        all_views = (
            scatterplot_views
            + ([obs_sets_view] if obs_sets_view is not None else [])
            + ([feature_list_view] if feature_list_view is not None else [])
        )
        rows_of_views = []

        for view_i, view in enumerate(all_views):
            if view_i % ncols == 0:
                # Start a new row.
                rows_of_views.append([view])
            else:
                rows_of_views[-1].append(view)
        
        # Combine rows.
        row_layouts = [hconcat(*row) for row in rows_of_views]
        vc.layout(vconcat(*row_layouts))

    else:
        # Single color or color=None: single scatterplot view.
        scatterplot_props = {} if title is None else { 'title': title }
        scatterplot_view = vc.add_view(vt.SCATTERPLOT, dataset=dataset).set_props(**scatterplot_props)
        control_views = []
        if has_color_in_obs:
            obs_sets_view = vc.add_view(vt.OBS_SETS, dataset=dataset)
            control_views.append(obs_sets_view)
        if has_color_in_var:
            feature_list_view = vc.add_view(vt.FEATURE_LIST, dataset=dataset)
            control_views.append(feature_list_view)
        
        # Link views.
        vc.link_views([scatterplot_view],
            ["embeddingType"],
            [basis_name]
        )

        if size is not None:
            vc.link_views([scatterplot_view],
                ["embeddingObsRadiusMode", "embeddingObsRadius"],
                ["manual", size]
            )
        
        if vmin is not None or vmax is not None:
            vc.link_views([scatterplot_view],
                ["featureValueColormapRange"],
                [
                    [
                        vmin if vmin is not None else 0.0,
                        vmax if vmax is not None else 1.0,
                    ]
                ]
            )
        
        if has_color_in_obs:
            vc.link_views([scatterplot_view, *control_views],
                ["obsColorEncoding"],
                ["cellSetSelection"]
            )
            vc.link_views([obs_sets_view],
                ["obsSetExpansion"],
                [[[obs_set_names[0]]]]
            )

        if has_color_in_var:
            vc.link_views([scatterplot_view, *control_views],
                ["featureSelection", "obsColorEncoding"],
                [[color[0]], "geneSelection"]
            )
        
        if color_map is not None and color_map in ["viridis", "plasma", "jet", "greys"]:
            vc.link_views([scatterplot_view],
                ["featureValueColormap"],
                [color_map]
            )

        # Layout views.
        if len(control_views) > 0:
            vc.layout(hconcat(scatterplot_view, vconcat(*control_views), split=[2, 1]))
        else:
            vc.layout(scatterplot_view)

    vw = _to_widget(vc)
    return vw


def spatial(adata, **kwargs):
    """
    This plotting function is deprecated since Scanpy version 1.11.0.

    :param AnnData adata: AnnData object.
    :param str color: Gene.
    :param str color_map: Color map (viridis, plasma, jet). Defaults to viridis.
    :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
    """

    warnings.warn("This plotting function is deprecated since Scanpy version 1.11.0.", DeprecationWarning)

    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sample_id = (list(adata.uns["spatial"].keys()))[0]
    
    color = kwargs.get("color", "")
    color_map = kwargs.get("color_map", "viridis")

    ncols = kwargs.get("ncols", 1)
    if ncols > 3:
        warnings.warn("To prevent plots from being too small, ncols should be ≤ 3.", UserWarning)
    
    
    output_img = join("data", "spatial.ome.zarr")
    output_adata = join("data", "spatial.anndata.zarr")

    x = 1/adata.uns['spatial'][sample_id]["scalefactors"]["tissue_hires_scalef"]
    y = x

        # Write img_arr to OME-Zarr.
    # Need to convert images from interleaved to non-interleaved (color axis should be first).
    img_hires = adata.uns['spatial'][sample_id]['images']['hires']
    img_arr = np.transpose(img_hires, (2, 0, 1))
    # Convert values from [0, 1] to [0, 255].
    img_arr *= 255.0
    
    # First, save the image to an OME-Zarr image format
    rgb_img_to_ome_zarr(img_arr, output_img, axes="cyx", chunks=(1, 256, 256), img_name="Image")
    # Second, save the AnnData object to Zarr format
    adata.obsm["spatial"] = adata.obsm["spatial"].astype("int32")

    adata_wrapper_dict = {
        "adata_path":output_adata,
        # obs_feature_matrix_path ?
        "coordination_values":{
            "obsType": 'cell',
            # 'featureType': 'gene'
            # "featureType": 'qualityMetric',
            # "featureValueType": 'value',
             # "featureValueType": 'exression'
            # obsLabelsType = null?
        }
    }

    if (type(color) == list and color[0] in adata.var.index) or (type(color) == str and color in adata.var.index): # gene
        # genes = kwargs["color"]
        # adata.var["genes"] = list(adata.var.index)
        # adata.var["in_color"] = adata.var["genes"].apply(lambda gene: True if gene in color else False)
        
        path = {"obs_feature_matrix_path": "X"}
        #path = {"feature_filter_path": ["var/in_color"]}
        new_coord_vals = {"featureType": 'gene', "featureValueType": 'expression'}
        # new_coord_vals = {"featureType": 'qualityMetric', "featureValueType": 'value'}
        adata_wrapper_dict.update(path)
        adata_wrapper_dict["coordination_values"].update(new_coord_vals)
        print(adata_wrapper_dict)
        
    elif (type(color) == list and color[0] in adata.obs.columns) or (type(color) == str and color in adata.obs.columns): # categorical
        adata.obs[color] = adata.obs[color].astype("float32")
        
        path = {"obs_feature_column_paths":[f"obs/{color}"]}
        new_coord_vals = {"obsType": 'cell', "featureType": 'qualityMetric', "featureValueType": 'value'}
        adata_wrapper_dict.update(path)
        adata_wrapper_dict["coordination_values"].update(new_coord_vals)
        print(adata_wrapper_dict)
    
    adata.write_zarr(output_adata, chunks=[adata.shape[0], VAR_CHUNK_SIZE])

    # obs_feature_column_paths=[f"obs/{color}"],
    # feature_filter_path=[f"obs/{color}"],
        
        
    vc = VitessceConfig(schema_version="1.0.18", name="sc.pl.spatial (deprecated)")
    dataset = vc.add_dataset(name='spatial data').add_object(
        AnnDataWrapper(
            adata_path=output_adata,
            obs_spots_path = "obsm/spatial", 
            obs_feature_matrix_path = "X"
            # obs_set_paths = ["obs/log1p_n_genes_by_counts"]
        )
    ).add_object(
        ImageOmeZarrWrapper(
            img_path=output_img,
            coordinate_transformations = [
                {
                    "type": 'translation',
                    "translation": [0, 0, 0],
                },
                {
                    "type": 'scale',
                    "scale": [1, x, y],
                    # [color, x, y]
                },
            ]
        )
    ).add_object(
        AnnDataWrapper(**adata_wrapper_dict)
    )
        # adata_path=output_adata,
        # obs_feature_column_paths=[f"obs/{color}"], # for numerical data
        # # obs_feature_matrix_path = []
        # coordination_values={
        #     "obsType": 'cell',
        #     "featureType": 'qualityMetric',
        #     "featureValueType": 'value',
        # }
    link_views_dict = {
        "obsType": 'cell',
        "featureSelection": [color],
        "obsColorEncoding": "geneSelection" # ??
    }

    if color in adata.var.index: # gene
        genes = vc.add_view(vt.FEATURE_LIST, dataset=dataset) #assumes featureType = gene
        
    if color in adata.obs.columns:
        histogram = vc.add_view(vt.FEATURE_VALUE_HISTOGRAM, dataset=dataset)

    # if type(color) == list and len(color) > 1:
    #     for i in range(0, len(color), 2):
    #         adata.obs[color[i]] = adata.obs[color[i]].astype("float32")
    #         path = {"obs_feature_column_paths":[f"obs/{color[i]}"]}
    #         new_coord_vals = {"obsType": 'cell', "featureType": 'qualityMetric', "featureValueType": 'value'}
    #         adata_wrapper_dict.update(path)
    #         adata_wrapper_dict["coordination_values"].update(new_coord_vals)
    #         # print(adata_wrapper_dict)
    #         spotLayer = {
    #         "spotLayer": CL([
    #             {
    #                 "obsType": "cell",
    #                 "spatialSpotRadius": 45, #might have to depend on scale factor as well
    #                 "featureValueColormap": color_map
    #             },
    #         ]) }
    #         spatial_view_1 = vc.add_view("spatialBeta", dataset=dataset)
    #         spatial_view_2 = vc.add_view("spatialBeta", dataset=dataset)
    
    #         vc.link_views_by_dict([spatial_view_1], **spotLayer, **link_views_dict)
    #         vc.link_views_by_dict([spatial_view_2],  **spotLayer, **link_views_dict)
            
    #         vc.layout(spatial_view_1 | spatial_view_2)
    
    #     vw = _to_widget(vc)
    #     return vw

    
    spatial_view = vc.add_view("spatialBeta", dataset=dataset)
    lc_view = vc.add_view("layerControllerBeta", dataset=dataset)

    
    if color in adata.obs.columns: # categorical
        new_vals = {"featureType": 'qualityMetric', "featureValueType": 'value'}
        link_views_dict.update(new_vals)
        print(link_views_dict)

    elif color in adata.var.index: # gene
        new_vals = {"featureType": 'gene', "featureValueType": 'expression'}
        link_views_dict.update(new_vals)
        print(link_views_dict)
        
    link_views_dict_without_feature_selection = {k:v for k, v in link_views_dict.items() if k != "featureSelection"}
    
    vc.link_views_by_dict([spatial_view, lc_view, genes],  {
        "spotLayer": CL([
            {
                "obsType": "cell",
                "spatialSpotRadius": 45, #might have to depend on scale factor as well
                "featureValueColormap": color_map
            },
        ]),
        **link_views_dict_without_feature_selection
    })

    vc.link_views([spatial_view, lc_view, genes], ["featureSelection"], [link_views_dict["featureSelection"]])
    
    vc.layout(spatial_view | (lc_view / (histogram if color in adata.obs.columns else genes)))
    
    vw = _to_widget(vc)
    return vw


# References:
# - https://github.com/scverse/scanpy/blob/cf8b46dea735c35a629abfaa2e1bab9047281e34/src/scanpy/plotting/_anndata.py#L1107
# - https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.heatmap.html
def heatmap(
        adata,
        var_names,
        groupby,
        *,
        gene_symbols=None,
        layer=None,
        swap_axes=False,
        vmin=None,
        vmax=None,
        **kwargs
    ):
    """
    Creates interactive heatmap.

    :param AnnData adata: AnnData object.
    :param list[str] var_names: List of genes.
    :param str groupby: Category group, a column in adata.obs.
    :param str gene_symbols: Key in adata.var to use for gene symbols, if different from adata.var.index.
    :param str layer: Layer in AnnData to use for expression values. Defaults to None, which uses adata.X.
    :param bool swap_axes: If True, transpose the heatmap. Defaults to False.
    :param float vmin: Minimum value for color map scaling.
    :param float vmax: Maximum value for color map scaling.
    :param str color_map: Color map (viridis, plasma, jet). Defaults to viridis.
    :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
    """

    # TODO: throw errors/warnings for unsupported parameters?

    vc =  VitessceConfig(schema_version="1.0.18", name='sc.pl.heatmap')

    # We need to make a copy since we will be modifying adata.
    # TODO: this will mean the object will always be re-written to disk
    # (and redundantly if the function is called multiple times),
    # and breaking the functionality of register_data_path() if the user had specified an existing on-disk location.
    # Workarounds for future work (TODO: create github issues):
    # - Create a new AnnData object (or CSV file) to store var/__ev_initial_feature_filter__ column only.
    # - In JS, implement a way to specify a list of features directly rather than via a path in the AnnData object.
    adata = adata.copy()

    # We need to create a boolean mask in adata.var to indicate which genes to show initially.
    # We name it using double underscores to reduce the chance of name collisions.
    adata.var["__ev_initial_feature_filter__"] = adata.var_names.to_series().apply(lambda gene: True if gene in var_names else False)

    group_name = groupby.capitalize()

    wrapper_params = {
        "obs_set_paths": [f"obs/{groupby}"],
        "obs_set_names": [group_name],
        "obs_feature_matrix_path": "X" if layer is None else f"layers/{layer}",
        "initial_feature_filter_path": "var/__ev_initial_feature_filter__",
    }

    if gene_symbols is not None:
        wrapper_params["feature_labels_path"] = f"var/{gene_symbols}"

    dataset = vc.add_dataset(name='heatmap data').add_object(AnnDataWrapper(
        **_get_adata_wrapper_params(adata),
        **wrapper_params,
    ))

    obs_sets_view = vc.add_view(vt.OBS_SETS, dataset=dataset)
    heatmap_view = vc.add_view(vt.HEATMAP, dataset=dataset).set_props(transpose=(not swap_axes))

    if "color_map" in kwargs:
        # TODO: color_map is not listed as a parameter for sc.pl.heatmap in Scanpy docs.
        # Should we still keep it here?
        color_map = kwargs["color_map"]
        vc.link_views([heatmap_view, obs_sets_view], 
            ["featureValueColormap"],
            [color_map]
        )
    
    vc.link_views([obs_sets_view], 
        ['obsSetExpansion'],
        [[[group_name]]]
    )

    if vmin is not None or vmax is not None:
        vc.link_views(
            [heatmap_view, obs_sets_view], 
            ["featureValueColormapRange"],
            [
                [
                    vmin if vmin is not None else 0.0,
                    vmax if vmax is not None else 1.0,
                ]
            ]
        )

    vc.layout(hconcat(heatmap_view, obs_sets_view, split=[2, 1]))
    vw = _to_widget(vc)
    return vw

# References:
# - https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.violin.html
# - https://github.com/scverse/scanpy/blob/cf8b46dea735c35a629abfaa2e1bab9047281e34/src/scanpy/plotting/_anndata.py#L735
def violin(
        adata,
        keys,
        groupby,
        *,
        log=False,
        stripplot=True,
        jitter=True,
        layer=None,
        order=None,
        **kwargs
    ):
    """
    Creates interactive violin plot.

    :param Anndata adata: AnnData object.
    :param list[str] keys: Gene IDs from adata.var.index.
    :param str groupby: Category group, a column in adata.obs.
    :param bool log: If True, apply log1p transformation to expression values before plotting. Defaults to False.
    :param bool stripplot: If True, add a strip plot alongside the violin plot. Defaults to True.
    :param bool jitter: If True, add jitter to the strip plot. Defaults to True.
    :param str layer: Layer in AnnData to use for expression values. Defaults to None, which uses adata.X.
    :param list[str] order: Order of categories in groupby.
    :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
    """
    vc =  VitessceConfig(schema_version="1.0.18", name='sc.pl.violin')
    
    if type(keys) == str:
        markers = [keys]
    elif type(keys) == list: 
        markers = keys
    else:
        raise ValueError("Parameter 'keys' must be a string or a list of strings.")

    group_name = groupby.capitalize()
    
    obs_set_selection = None
    if order is not None:
        obs_set_selection = [
            # Construct the obsSetPaths in the order specified by `order`.
            [group_name, item] for item in order
        ]
    

    dataset = vc.add_dataset(name='violin data').add_object(AnnDataWrapper(
        **_get_adata_wrapper_params(adata),
        obs_set_paths=[f"obs/{groupby}"],
        obs_set_names=[group_name],
        obs_feature_matrix_path="X" if layer is None else f"layers/{layer}",
    ))

    # Case: multiple genes. Create one violin plot view per gene.
    if len(markers) > 1:
        feature_list_view = vc.add_view(vt.FEATURE_LIST, dataset=dataset)
        obs_sets_view = vc.add_view(vt.OBS_SETS, dataset=dataset)
        # Create one violin plot per gene, and add to the list violin_views.
        violin_views = []
        for gene_i, gene in enumerate(markers):
            violin_view = vc.add_view('obsSetFeatureValueDistribution', dataset=dataset).set_props(jitter=(stripplot and jitter))

            if gene_i == 0:
                # We want to link the first violin plot to the feature list, but not the others.
                vc.link_views([violin_view, feature_list_view], 
                    ["featureSelection"],
                    [[gene]]
                )
            else:
                vc.link_views([violin_view], 
                    ["featureSelection"],
                    [[gene]]
                )
            
            vc.link_views([violin_view, obs_sets_view], 
                ['obsSetExpansion'],
                [[[group_name]]]
            )

            if obs_set_selection is not None:
                vc.link_views([violin_view, obs_sets_view],
                    ["obsSetSelection"],
                    [obs_set_selection]
                )

            if log:
                vc.link_views([violin_view],
                    ["featureValueTransform"],
                    ["log1p"]
                )
            violin_views.append(violin_view)

        # Layout all of the views.
        vc.layout(hconcat(vconcat(*violin_views), vconcat(feature_list_view, obs_sets_view), split = [2, 1]))
    else:
        # Case: single gene. Create one violin plot view.
        feature_list_view = vc.add_view(vt.FEATURE_LIST, dataset=dataset)
        obs_sets_view = vc.add_view(vt.OBS_SETS, dataset=dataset)
        violin_view = vc.add_view('obsSetFeatureValueDistribution', dataset=dataset).set_props(jitter=(stripplot and jitter))

        vc.link_views([violin_view, feature_list_view], 
            ["featureSelection"],
            [markers]
        )
        vc.link_views([violin_view, obs_sets_view], 
            ['obsSetExpansion'],
            [[[group_name]]]
        )

        if obs_set_selection is not None:
            vc.link_views([violin_view, obs_sets_view],
                ["obsSetSelection"],
                [obs_set_selection]
            )

        if log:
            vc.link_views([violin_view],
                ["featureValueTransform"],
                ["log1p"]
            )

        vc.layout(hconcat(violin_view, vconcat(feature_list_view, obs_sets_view), split=[2, 1]))

    vw = _to_widget(vc)
    return vw

# References:
# - https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.dotplot.html
# - https://github.com/scverse/scanpy/blob/cf8b46dea735c35a629abfaa2e1bab9047281e34/src/scanpy/plotting/_dotplot.py#L844
def dotplot(
        adata,
        var_names,
        groupby,
        *,
        expression_cutoff=None,
        title=None,
        gene_symbols=None,
        layer=None,
        cmap=None,
        swap_axes=False,
        **kwargs,
    ):
    """
    Creates interactive dotplot.

    :param AnnData adata: AnnData object.
    :param list[str] var_names: List of gene IDs from adata.var.index.
    :param str groupby: Category group, a column in adata.obs.
    :param float expression_cutoff: Expression cutoff for dot display.
    :param str title: Title of the plot.
    :param str gene_symbols: Key in adata.var to use for gene symbols, if different from adata.var.index.
    :param str layer: Layer in AnnData to use for expression values. Defaults to None, which uses adata.X.
    :param str cmap: Color map (viridis, plasma, jet). Defaults to viridis.
    :param bool swap_axes: If True, transpose the dotplot. Defaults to False.
    :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
    """

    vc = VitessceConfig(schema_version="1.0.18", name='sc.pl.dotplot')

    group_name = groupby.capitalize()

    wrapper_params = {
        "obs_set_paths": [f"obs/{groupby}"],
        "obs_set_names": [group_name],
        "obs_feature_matrix_path": "X" if layer is None else f"layers/{layer}",
    }

    if gene_symbols is not None:
        wrapper_params["feature_labels_path"] = f"var/{gene_symbols}"
    
    dataset = vc.add_dataset('dotplot data').add_object(AnnDataWrapper(
        **_get_adata_wrapper_params(adata),
        **wrapper_params,
    )).add_object(AnnDataWrapper(
        **_get_adata_wrapper_params(adata),
    ))

    obs_sets_view = vc.add_view(vt.OBS_SETS, dataset=dataset)
    feature_list_view = vc.add_view(vt.FEATURE_LIST, dataset=dataset).set_props(enableMultiSelect=True)
    dot_plot_view = vc.add_view('dotPlot', dataset=dataset).set_props(title=title, transpose=swap_axes)

    vc.link_views([obs_sets_view], 
        ['obsSetExpansion'],
        [[[group_name]]]
    )
    
    if var_names is not None:
        vc.link_views([dot_plot_view, feature_list_view],
            ["featureSelection"],
            [var_names]
        )
    
    if expression_cutoff is not None:
        vc.link_views([dot_plot_view],
            ["featureValuePositivityThreshold"],
            [expression_cutoff]
        )
    
    if cmap is not None and cmap in ["viridis", "plasma", "jet", "greys"]:
        # TODO: support a "reds" colormap to reflect the scanpy default.
        vc.link_views([dot_plot_view],
            ["featureValueColormap"],
            [cmap]
        )

    vc.layout(hconcat(dot_plot_view, vconcat(feature_list_view, obs_sets_view), split=[2, 1]))
    vw = _to_widget(vc)
    return vw

def _monkeypatch(cls, func):
    """
    Modifies behavior of the class to replace a function.

    :param any cls: Class to be modified. Expected to be sc.pl class.
    :param any func: function to be replaced. Expected to be plotting function from sc.pl class.
    """

    func_name = func.__name__
    orig_func_name = f"_orig_{func_name}"
    if not hasattr(cls, orig_func_name):
        orig_func = getattr(cls, func_name)
        setattr(cls, orig_func_name, orig_func)
    setattr(cls, func_name, func)

def _undo_monkeypatch(cls, func_name):
    """
    Restores the original behavior of the class.
    """
    orig_func_name = f"_orig_{func_name}"
    if hasattr(cls, orig_func_name):
        orig_func = getattr(cls, orig_func_name)
        setattr(cls, func_name, orig_func)

def _monkeypatch_spatialdata():
    """
    Replaces behavior of SpatialData.pl class with VitesscePlotAccessor.
    """
    VitesscePlotAccessor._is_enabled = True

    if not hasattr(SpatialData, "pl"):
        raise ValueError("The accessor SpatialData.pl does not yet exist. Please import spatialdata_plot first.")
    if not hasattr(SpatialData, '_orig_pl'):
        # Not yet monkeypatched.
        setattr(SpatialData, '_orig_pl', _CachedAccessor('_orig_pl', SpatialData.pl))
        setattr(SpatialData, 'pl', _CachedAccessor('pl', VitesscePlotAccessor))
    else:
        print("Warning: SpatialData.pl has already been monkeypatched.")
    
def _undo_monkeypatch_spatialdata():
    """
    Restores the original behavior of SpatialData.pl.
    """
    VitesscePlotAccessor._is_enabled = False

    if not hasattr(SpatialData, "pl"):
        raise ValueError("The accessor SpatialData.pl does not yet exist. Please import spatialdata_plot first.")

    if hasattr(SpatialData, '_orig_pl'):
        # Has already been monkeypatched. Undo.
        setattr(SpatialData, 'pl', _CachedAccessor('pl', SpatialData._orig_pl))
        delattr(SpatialData, '_orig_pl')
    else:
        print("Warning: SpatialData.pl has not been monkeypatched yet.")

def configure_plots(disable_plots=None, enable_plots=None): 
    """
    Deactivates and reactivates interactive Vitessce plots.

    :param list[str] disable_plots: List of plots.
    :param list[str] enable_plots: List of plots.
    """

    SCANPY_PLOTTING_FUNCTIONS = {
        "embedding": embedding,
        "umap": umap,
        "pca": pca,
        "tsne": tsne,
        "diffmap": diffmap,
        "spatial": spatial,
        "dotplot": dotplot,
        "heatmap": heatmap,
        "violin": violin
    }
    ALL_PLOT_NAMES = list(SCANPY_PLOTTING_FUNCTIONS.keys()) + ["spatialdata-plot"]

    if type(enable_plots) == str or type(disable_plots) == str:
        raise RuntimeError("Expected enable_plots/disable_plots to be a list of string.")

    if type(enable_plots) == list and type(disable_plots) == list:
        if any(plot in enable_plots for plot in disable_plots):
            raise RuntimeError("Plots cannot be in enable_plots and disable_plots simultaneously.")
    
    # By default, enable all plots.
    if disable_plots is None:
        disable_plots = []

    if enable_plots is None:
        enable_plots = list(set(ALL_PLOT_NAMES) - set(disable_plots))

    for plot, func in SCANPY_PLOTTING_FUNCTIONS.items():
        if plot in enable_plots:
            _monkeypatch(sc.pl, func)
        elif plot in disable_plots:
            _undo_monkeypatch(sc.pl, plot)
            print(f"Deactivated Vitessce {plot}")

    if "spatialdata-plot" in enable_plots:
        _monkeypatch_spatialdata()
    elif "spatialdata-plot" in disable_plots:
        _undo_monkeypatch_spatialdata()
        print("Deactivated Vitessce spatialdata-plot")


# Convenience functions for enabling/disabling.
ALL_PLOTS = [
    "embedding",
    "umap",
    "pca",
    "tsne",
    "diffmap",
    "spatial",
    "dotplot",
    "heatmap",
    "violin",
    "spatialdata-plot"
]

def enable_plots(plots=None):
    """
    Activates interactive Vitessce plots.

    :param list[str] plots: List of plots to enable. If None, enables all plots.
    """

    if plots is None:
        plots = ALL_PLOTS
    
    configure_plots(enable_plots=plots)

def disable_plots(plots=None):
    """
    Deactivates interactive Vitessce plots.

    :param list[str] plots: List of plots to disable. If None, disables all plots.
    """

    if plots is None:
        plots = ALL_PLOTS

    configure_plots(disable_plots=plots)