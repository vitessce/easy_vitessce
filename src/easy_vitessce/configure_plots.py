from vitessce import (
    VitessceConfig,
    AnnDataWrapper,
    SpatialDataWrapper,
    ImageOmeZarrWrapper,
    CoordinationLevel as CL,
    Component as cm,
    vconcat,
    hconcat
)

import os
import shutil
from os.path import join
import scanpy as sc

from vitessce.data_utils import (
    VAR_CHUNK_SIZE,
    rgb_img_to_ome_zarr
)

import pandas as pd
import numpy as np
import warnings

from anndata import (AnnData, read_h5ad)

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
  :param str color: Gene or category group.
  :param str color_map: Color map (viridis, plasma, jet).
  :param float or int size: Size of dots.
  :param bool include_gene_list: If a list of genes is passed in, True will add a gene list for the last plot. False by default.
  :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
  """
  return embedding(adata, basis="umap", **kwargs)

def tsne(adata, **kwargs):
  """
  Creates interactive t-SNE plot.

  :param AnnData adata: AnnData object.
  :param str color: Gene or category group.
  :param str color_map: Color map (viridis, plasma, jet).
  :param (float or int) size: Size of dots.
  :param bool include_gene_list: If a list of genes is passed in, True will add a gene list for the last plot. False by default.
  :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
  """
  return embedding(adata, basis="tsne", **kwargs)

def pca(adata, **kwargs):
  """
  Creates interactive PCA plot.

  :param AnnData adata: AnnData object.
  :param str color: Gene or category group.
  :param str color_map: Color map (viridis, plasma, jet).
  :param (float or int) size: Size of dots.
  :param bool include_gene_list: If a list of genes is passed in, True will add a gene list for the last plot. False by default.
  :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
  """
  return embedding(adata, basis="pca", **kwargs)

def embedding(adata, basis, **kwargs):
    """
    Creates interactive versions of UMAP, PCA, t-SNE plots.

    :param AnnData adata: AnnData object.
    :param str basis: Name of plot (umap, pca, or tsne).
    :param str color: Gene or categorical label.
    :param str color_map: Color map (viridis, plasma, jet). Defaults to viridis.
    :param (float or int) size: Size of dots.
    :param bool include_gene_list: If a list of genes is passed in, True will add a gene list for the last plot. False by default.
    :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 

    """
    basis = basis
    basis_name = "t-SNE" if basis == "tsne" else basis.upper()
    basis_obsm_key = f"X_{basis}"
    adata = adata

    if basis_obsm_key not in adata.obsm:
        raise ValueError(f"{basis_name} coordinates not found in adata.obsm.")
    
    ncols = kwargs.get("ncols", 1)
    if ncols > 3:
        warnings.warn("To prevent plots from being too small, ncols should be ≤ 3.")

    include_genes = kwargs.get("include_gene_list", False) # Note: Not a parameter in Scanpy.
    include_cells = kwargs.get("include_cell_sets", False) # Note: Not a parameter in Scanpy.
    
    color = kwargs.get("color", "")
    color_map = kwargs.get("color_map", "viridis")
    size  = kwargs.get("size", 2.5)

    if "color_map" in kwargs.keys():
        if "plasma" not in kwargs["color_map"].lower() and "viridis" not in kwargs["color_map"].lower() and "jet" not in kwargs["color_map"].lower():
            print("Invalid color_map. Supported color_maps: plasma, viridis, jet.")
            color_map = "viridis"
            
    coordination_types = ["embeddingObsRadiusMode", "embeddingObsRadius", "embeddingObsOpacityMode", "obsColorEncoding"]
    coordination_values = ["manual", size, "manual"]
    
    if (type(color) == str and len(color) != 0) or len(color) == 1:
        # color = kwargs.get("color") if type(kwargs.get("color")) == str else kwargs.get("color")[0]
        
        if color in adata.obs.columns:
            # `color` is the name of a categorical `obs` dataframe column
            new_coord_values = ["cellSetSelection", [color if type(color) == str else color[0]], color_map]
            new_coord_types = ["featureSelection", "featureValueColormap"]
            coordination_values.extend(new_coord_values)
            coordination_types.extend(new_coord_types)
            
            
        elif color in adata.var.index:
            # `color` is a gene name
            new_coord_values = ["geneSelection", [color if type(color) == str else color[0]], color_map]
            new_coord_types = ["featureSelection", "featureValueColormap"]
            coordination_values.extend(new_coord_values)
            coordination_types.extend(new_coord_types) 
        else:
          raise ValueError("color not found in AnnData object")
        
    elif type(kwargs.get("color")) == list and len(kwargs.get("color")) > 1: 
        color = kwargs.get("color", [])
        
    vc =  VitessceConfig(schema_version="1.0.15", name=basis_name)

    adata_wrapper_dict = {
        **_get_adata_wrapper_params(adata),
        'obs_embedding_paths':[f"obsm/{basis_obsm_key}"],
        'obs_embedding_names':[basis_name],
        'obs_feature_matrix_path':"X"
    }

    if len(color) == 0:
        dataset = vc.add_dataset(name='data').add_object(AnnDataWrapper(**adata_wrapper_dict))
        coordination_types.remove("obsColorEncoding")
        coordination_types.append("featureValueColormap")
        coordination_values.append(color_map)
        mapping = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping=basis_name) # mapping value corresponds to one of the obs_embedding_names values.
        view_list = vc.add_view(cm.FEATURE_LIST, dataset=dataset) # default to feature_list if no color is provided

        vc.link_views(
            [mapping, view_list], 
            coordination_types, # https://vitessce.io/docs/coordination-types/
            coordination_values
        )

        vc.layout(mapping | view_list)
       
        vw = vc.widget()
        return vw
    
    elif type(color) == list and len(color) > 1: 
        
        if color[0] in adata.var.index: # gene names
            dataset = vc.add_dataset(name='data').add_object(AnnDataWrapper(**adata_wrapper_dict))
        
            mapping_dict = {}
            for i in range(0, len(color), ncols):
                row = color[i:i + ncols] # not inclusive
                for gene in row:
                    mapping_dict[f"{gene}"] = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping=basis_name).set_props(title=f"{basis_name} ({gene})")
            
            scatterplots = list(mapping_dict.values())
            vc.link_views(scatterplots,
                ["embeddingZoom", "embeddingTargetX", "embeddingTargetY"],
                [None, None, None]
            )
                                                                
            for key, value in mapping_dict.items():
               vc.link_views(
                    [value], 
                    ["featureSelection", "obsColorEncoding", "embeddingObsRadiusMode", "embeddingObsRadius", "featureValueColormap"], # https://vitessce.io/docs/coordination-types/
                    [[key], "geneSelection", "manual", size, color_map]
                )
            
            
            values = [mapping_dict[key] for key in mapping_dict.keys()]
    
            cols = []
            for i in range(0, len(values), ncols):
                cols.append(values[i:i+ncols])
            #print(cols)
    
            if len(values) % ncols != 0:
                last_row = [values[len(values) - len(values)%ncols]]
                last_key = (list(mapping_dict.keys()))[-1]
    
                if include_genes or include_cells: # hmmm
                    last_gene_list = vc.add_view(cm.FEATURE_LIST, dataset=dataset)
                    last_cell_list = vc.add_view(cm.OBS_SETS, dataset=dataset)
                    
                    cols[len(cols)-1].append(last_gene_list)
                    
                    vc.link_views(
                        [last_row[0], last_gene_list], 
                        ["featureSelection", "obsColorEncoding", "embeddingObsRadiusMode", "embeddingObsRadius", "featureValueColormap"], # https://vitessce.io/docs/coordination-types/
                        [[last_key], "geneSelection", "manual", size, color_map]
                    )

                    #last_row.append(last_gene_list)
                else:
                    vc.link_views(
                        [last_row[0]], 
                        ["featureSelection", "obsColorEncoding", "embeddingObsRadiusMode", "embeddingObsRadius", "featureValueColormap"], # https://vitessce.io/docs/coordination-types/
                        [[last_key], "geneSelection", "manual", size, color_map]
                    )
                
    
                vc.layout((vconcat(*[hconcat(*row) for row in cols])))
                    
                    
            else:
                last_mapping = (list(mapping_dict.values()))[-1]
                last_key = (list(mapping_dict.keys()))[-1]
    
                if include_genes:
                    genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset)
                    cols.append([genes])
                    vc.link_views(
                        [last_mapping, genes], 
                        ["featureSelection", "obsColorEncoding", "embeddingObsRadiusMode", "embeddingObsRadius", "featureValueColormap"], # https://vitessce.io/docs/coordination-types/
                        [[last_key], "geneSelection", "manual", size, color_map]
                    )
                    vc.layout(genes)
    
                vc.layout((vconcat(*[hconcat(*row) for row in cols])))
                
        elif color[0] in adata.obs.columns: # categorical, layout shouldn't be side-by-side 
            for obs in color:
                # print(obs)
                new_obs_paths_names = {'obs_set_paths': [f"obs/{obs}"], 'obs_set_names':[obs]} # is the problem here?
                adata_wrapper_dict.update(new_obs_paths_names)
                # print(adata_wrapper_dict)
                
                dataset = vc.add_dataset(name='tsne data').add_object(AnnDataWrapper(**adata_wrapper_dict))
                
                scatterplot = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping=basis_name)
                obs_view = vc.add_view(cm.OBS_SETS, dataset=dataset)
            
                vc.link_views(
                    [scatterplot, obs_view], 
                    ["obsColorEncoding", "embeddingObsRadiusMode", "embeddingObsRadius", "embeddingObsOpacityMode", "obsSetSelection", "obsSetColor", 'obsSetExpansion'],
                    ["cellSetSelection", "manual", size, "manual", None, None, [[obs]]]
                )
            
                vc.layout(scatterplot | obs_view)
                
            
        vw = _to_widget(vc)
        return vw

    else: # one color
        if color in adata.var.index:
            dataset = vc.add_dataset(name='data').add_object(AnnDataWrapper(**adata_wrapper_dict))
        
        elif color in adata.obs.columns:
            first_color = color if type(color) == str else color[0]
            obs_paths_names = {'obs_set_paths':[f"obs/{first_color}"], 'obs_set_names':[first_color.capitalize()]}
            adata_wrapper_dict.update(obs_paths_names)
            
        dataset = vc.add_dataset(name='data').add_object(AnnDataWrapper(**adata_wrapper_dict))
        
        mapping = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping=basis_name) # mapping value corresponds to one of the obs_embedding_names values.
        view_list = vc.add_view(cm.FEATURE_LIST if color in adata.var.index else cm.OBS_SETS, dataset=dataset) #change dimensions?
    

        vc.link_views(
            [mapping, view_list], 
            coordination_types, # https://vitessce.io/docs/coordination-types/
            coordination_values
        )

        vc.layout(mapping | view_list)
       
    vw = _to_widget(vc)
    return vw

def spatial(adata, **kwargs):
    """
    Creates interactive spatial plot. Similar syntax to Scanpy's spatial plot.

    :param AnnData adata: AnnData object.
    :param str color: Gene.
    :param str color_map: Color map (viridis, plasma, jet). Defaults to viridis.
    :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
    """
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sample_id = (list(adata.uns["spatial"].keys()))[0]
    
    color = kwargs.get("color", "")
    color_map = kwargs.get("color_map", "viridis")

    ncols = kwargs.get("ncols", 1)
    if ncols > 3:
        warnings.warn("To prevent plots from being too small, ncols should be ≤ 3.")
    
    
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
        
        
    vc = VitessceConfig(schema_version="1.0.17", name="AnnData with image")
    dataset = vc.add_dataset("My dataset").add_object(
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
        genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset) #assumes featureType = gene
        
    if color in adata.obs.columns:
        histogram = vc.add_view(cm.FEATURE_VALUE_HISTOGRAM, dataset=dataset)

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

def heatmap(adata, **kwargs):
    """
    Creates interactive heatmap.

    :param AnnData adata: AnnData object.
    :param str groupby: Category group.
    :param list[str] var_names: List of genes.
    :param str color_map: Color map (viridis, plasma, jet). Defaults to viridis.
    :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
    """
    vc =  VitessceConfig(schema_version="1.0.15", name='heatmap')
    adata = adata
    color_map = kwargs.get("color_map", "viridis")

    if "groupby" in kwargs:
        groupby = kwargs["groupby"]

    if "var_names" in kwargs.keys():
        markers = kwargs["var_names"]
        adata.var["genes"] = list(adata.var.index)
        adata.var["in_markers"] = adata.var["genes"].apply(lambda gene: True if gene in markers else False)


    dataset = vc.add_dataset(name='data').add_object(AnnDataWrapper(
        **_get_adata_wrapper_params(adata),
        obs_set_paths=[f"obs/{groupby}"],
        obs_set_names=[groupby],
        # obs_embedding_paths=["obsm/X_umap"],
        # obs_embedding_names=["UMAP"],
        obs_feature_matrix_path="X",
        initial_feature_filter_path="var/in_markers"
    ))

    cells = vc.add_view(cm.OBS_SETS, dataset=dataset)
    heatmap = vc.add_view(cm.HEATMAP, dataset=dataset).set_props(transpose=True)

    vc.link_views(
        [heatmap, cells], 
        ["featureValueColormap"],
        [color_map]
    )
            
    vc.layout(heatmap | cells)

    vw = _to_widget(vc)
    return vw

def violin(adata, groupby,**kwargs):
    """
    Creates interactive violin plot.

    :param Anndata adata: AnnData object.
    :param str groupby: Category group.
    :param list[str] keys: Genes.
    :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
    """
    vc =  VitessceConfig(schema_version="1.0.15", name='heatmap')
    adata = adata
    groupby = groupby

    if "keys" not in kwargs.keys():
        markers = []
    
    if type(kwargs.get("keys")) == str:
        markers = [kwargs.get("keys", [])]
    elif type(kwargs.get("keys")) == list: 
        markers = kwargs.get("keys", [])
    
    dataset = vc.add_dataset(name='data').add_object(AnnDataWrapper(
        **_get_adata_wrapper_params(adata),
        obs_set_paths=[f"obs/{groupby}"],
        obs_set_names=[groupby],
        obs_feature_matrix_path="X"
    ))

    if type(markers) == list and len(markers) > 1:
        for gene in markers:
            genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset)
            cells = vc.add_view(cm.OBS_SETS, dataset=dataset)
            violin = vc.add_view('obsSetFeatureValueDistribution', dataset=dataset, uid=f'violin-plot-{gene}')
            vc.link_views(
                [violin, genes, cells], 
                ["featureSelection", "obsSetSelection", 'obsSetExpansion'],
                [[gene], None, [[gene]]]
            )
            vc.layout(hconcat(violin, genes, cells, split = [2,1,1]))
    else:
        genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset)
        cells = vc.add_view(cm.OBS_SETS, dataset=dataset)
        violin = vc.add_view('obsSetFeatureValueDistribution', dataset=dataset, uid='violin-plot')

        if "keys" in kwargs.keys():
            # print(markers)
            vc.link_views(
                [violin, genes, cells], 
                ["featureSelection", 'obsSetExpansion'],
                [markers, [markers] if type(markers) == list else [[markers]]]
            )
        
        vc.layout(violin | genes / cells)

    vw = _to_widget(vc)
    return vw

def dotplot(adata, groupby, **kwargs):
    """
    Creates interactive dotplot.

    :param AnnData adata: AnnData object.
    :param str groupby: Category group.
    :param list[str] var_names: List of genes.
    :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
    """
    markers = kwargs.get("var_names")

    vc = VitessceConfig(schema_version="1.0.17", name='dotplot data')
    
    dataset = vc.add_dataset('dotplot data').add_object(AnnDataWrapper(
        **_get_adata_wrapper_params(adata),
        obs_set_paths=[f"obs/{groupby}"],
        obs_set_names=[groupby],
        obs_feature_matrix_path="X",
    )).add_object(AnnDataWrapper(
        **_get_adata_wrapper_params(adata),
    ))


    obsSets = vc.add_view(cm.OBS_SETS, dataset=dataset)
    featureList = vc.add_view(cm.FEATURE_LIST, dataset=dataset).set_props(enableMultiSelect=True)
    dotPlot = vc.add_view('dotPlot', dataset=dataset)
    
    if markers is not None:
        vc.link_views(
            [dotPlot, featureList], 
            ["featureSelection"],
            [markers]
        )
    
    vc.layout(dotPlot | featureList / (obsSets))
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
        "spatial": spatial,
        "dotplot": dotplot,
        "heatmap": heatmap,
        "violin": violin
    }
    ALL_PLOT_NAMES = list(SCANPY_PLOTTING_FUNCTIONS.keys()) + ["spatialdata-plot"]

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
