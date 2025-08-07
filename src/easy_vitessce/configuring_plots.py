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

import spatialdata 
from spatialdata import SpatialData
from xarray.core.extensions import _CachedAccessor
from spatialdata_plot.pl.basic import PlotAccessor

from easy_vitessce.VitessceSpatialData import VitessceSpatialData

def _create_zarr_filepath(adata, plot_type):
    """
    Creates Zarr filepath for AnnData object. Prevents creating multiple files with the same name.

    :param AnnData adata: AnnData object.
    :param str plot_type: plot type.
    :returns: Zarr filepath.
    """ 
    zarr_filepath = join("data", f"{plot_type}_file.zarr")
    if os.path.exists(zarr_filepath) and os.path.isdir(zarr_filepath):
            shutil.rmtree(zarr_filepath)   
    adata.write_zarr(zarr_filepath, chunks=[adata.shape[0], VAR_CHUNK_SIZE])
    return zarr_filepath


def embedding(adata, basis, **kwargs):
    """
    Creates interactive versions of UMAP, PCA, t-SNE plots.

    :param AnnData adata: AnnData object.
    :param str basis: Name of plot (umap, pca, or tsne).
    :param str color: Gene.
    :param str color_map: Color map (viridis, plasma, jet). Defaults to viridis.
    :param (float or int) size: Size of dots.
    :param bool include_gene_list: If a list of genes is passed in, True will add a gene list for the last plot. False by default.
    :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 

    """
    basis = basis
    adata = adata

    ncols = kwargs.get("ncols", 1)
    if ncols > 3:
        warnings.warn("To prevent plots from being too small, ncols should be ≤ 3.")

    include_genes = kwargs.get("include_gene_list", False)
    
    if "color" not in kwargs.keys():
        color = ""

    if type(kwargs.get("color")) == str:
        color = kwargs.get("color", "")
    elif type(kwargs.get("color")) == list: 
        color = kwargs.get("color", [])
        
    color_map = kwargs.get("color_map", "viridis")
    size  = kwargs.get("size", 2.5)

    if "color_map" in kwargs.keys():
        if "plasma" not in kwargs["color_map"].lower() and "viridis" not in kwargs["color_map"].lower() and "jet" not in kwargs["color_map"].lower():
            print("Invalid color_map. Supported color_maps: plasma, viridis, jet. Default set to plasma.")
            color_map = "plasma"
            
    zarr_filepath = _create_zarr_filepath(adata, "embedding")
    
    
    if basis == "umap" or basis == "pca":
        vc =  VitessceConfig(schema_version="1.0.15", name='UMAP' if basis=="umap" else "PCA")
        
        dataset = vc.add_dataset(name='data').add_object(AnnDataWrapper(
            adata_store=zarr_filepath,
            obs_embedding_paths=["obsm/X_umap" if basis == "umap" else "obsm/X_pca"],
            obs_embedding_names=["UMAP" if basis == "umap" else "PCA"],
            obs_feature_matrix_path="X"
        ))

        if type(color) == list and len(color) > 1: 
            mapping_dict = {}
            #split = []
            #split.extend([1] * ncols)
            for i in range(0, len(color), ncols):
                row = color[i:i + ncols]
                #print(row) # not inclusive, used in layout later???
                for gene in row:
                    mapping_dict[f"{gene}"] = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="UMAP" if basis=="umap" else "PCA").set_props(title=f"{basis.upper()} ({gene})")
            
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
                # print(last_key)
                # print(f"last row: {last_row}")
                # print(type(last_row))
                #last_split = []
                #last_split.extend([1] * (len(last_row)+1)) # +1 accounts for gene list

                if include_genes:
                    last_gene_list = vc.add_view(cm.FEATURE_LIST, dataset=dataset)
                    
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

                #first_rows = values[0:len(values)- len(values)%ncols]

                #for mapping in mapping_dict.keys():
                    #vc.layout(vconcat((hconcat(*first_rows)), hconcat(*last_row)))
                    
                    
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
       

        else:
 
            mapping = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="UMAP" if basis=="umap" else "PCA") # mapping value corresponds to one of the obs_embedding_names values.
            genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset) #change dimensions?
        
        
            vc.link_views(
                [mapping, genes], 
                ["featureSelection", "obsColorEncoding", "embeddingObsRadiusMode", "embeddingObsRadius", "featureValueColormap"], # https://vitessce.io/docs/coordination-types/
                [[color[0] if type(color) == list else color], "geneSelection", "manual", size , color_map]
            )
            vc.layout(mapping | genes);
           
        vw = vc.widget()
        return vw


    elif basis == "tsne":
        if "X_tsne" not in adata.obsm:
            sc.tl.tsne(adata, random_state=1)
            zarr_filepath = _create_zarr_filepath(adata, "embedding")
       
        vc =  VitessceConfig(schema_version="1.0.15", name='t-SNE')
        
        if (type(color) == str) or (type(color) == list and len(color) == 1):
            dataset = vc.add_dataset(name='tsne data').add_object(AnnDataWrapper(
                    adata_path=zarr_filepath,
                    obs_set_paths=[f"obs/{color}" if type(color) == str else f"obs/{color[0]}"],
                    obs_set_names=["color"],
                    obs_embedding_paths=["obsm/X_tsne"],
                    obs_embedding_names=["t-SNE"],
                    obs_feature_matrix_path="X"
                ))
            
            tsne = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="t-SNE")
            cells = vc.add_view(cm.OBS_SETS, dataset=dataset)
        
            vc.link_views(
            [tsne, cells], 
            ["obsColorEncoding", "embeddingObsRadiusMode", "embeddingObsRadius", "embeddingObsOpacityMode", "embeddingObsOpacity"],
            ["cellSetSelection", "manual", size, "manual", 0.6]
        )
        
            vc.layout(tsne | cells)
            

        else:
            for obs in color:
                # print(obs)
                dataset = vc.add_dataset(name='tsne data').add_object(AnnDataWrapper(
                    adata_store=zarr_filepath,
                    obs_set_paths=[f"obs/{obs}"],
                    obs_set_names=[obs],
                    obs_embedding_paths=["obsm/X_tsne"],
                    obs_embedding_names=["t-SNE"],
                    obs_feature_matrix_path="X"
                ))
            
                tsne = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="t-SNE")
                obs = vc.add_view(cm.OBS_SETS, dataset=dataset)
            
                vc.link_views(
                [tsne, obs], 
                ["obsColorEncoding", "embeddingObsRadiusMode", "embeddingObsRadius", "embeddingObsOpacityMode", "embeddingObsOpacity"],
                ["cellSetSelection", "manual", size, "manual", 0.6]
            )
            
                vc.layout(tsne | obs)
                
            
        vw = vc.widget()
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
    adata.obs[color] = adata.obs[color].astype("float32")
    adata.write_zarr(output_adata, chunks=[adata.shape[0], VAR_CHUNK_SIZE])

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
    
        ).add_object(AnnDataWrapper(
        adata_path=output_adata,
        obs_feature_column_paths=[f"obs/{color}"],
        coordination_values={
            "obsType": 'cell',
            "featureType": 'qualityMetric',
            "featureValueType": 'value',
        }
    ))
    
    spatial_view = vc.add_view("spatialBeta", dataset=dataset)
    lc_view = vc.add_view("layerControllerBeta", dataset=dataset)
    # genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset) #assumes featureType = gene
    histogram = vc.add_view(cm.FEATURE_VALUE_HISTOGRAM, dataset=dataset)
    
    vc.link_views_by_dict([spatial_view, lc_view],  {
        "spotLayer": CL([
            {
                "obsType": "cell",
                "spatialSpotRadius": 45, #might have to depend on scale factor as well
                "featureValueColormap": color_map
            },
        ]),
            "obsType": 'cell',
        "featureType": 'qualityMetric',
        "featureValueType": 'value',
        "featureSelection": [color],
        "obsColorEncoding": "geneSelection"
    })
    
    vc.layout(spatial_view | lc_view / histogram)
    
    vw = vc.widget()
    return vw

def heatmap(adata, **kwargs):
    """
    Creates interactive heatmap.

    :param AnnData adata: AnnData object.
    :param str groupby: Category group.
    :param list[str] markers: List of genes.
    :param str color_map: Color map (viridis, plasma, jet). Defaults to viridis.
    :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
    """
    vc =  VitessceConfig(schema_version="1.0.15", name='heatmap')
    adata = adata
    color_map = kwargs.get("color_map", "viridis")

    if "groupby" in kwargs:
        groupby = kwargs["groupby"]

    if "markers" in kwargs:
        markers = kwargs["markers"]
        adata.var["genes"] = list(adata.var.index)
        adata.var["in_markers"] = adata.var["genes"].apply(lambda gene: True if gene in markers else False)

    zarr_filepath = _create_zarr_filepath(adata, "heatmap")

    dataset = vc.add_dataset(name='data').add_object(AnnDataWrapper(
            adata_path=zarr_filepath,
            obs_set_paths=[f"obs/{groupby}"],
            obs_set_names=["cell type"],
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

    vw = vc.widget()
    return vw

def violin(adata, groupby,**kwargs):
    """
    Creates interactive violin plot.

    :param Anndata adata: AnnData object.
    :param str groupby: Category group.
    :param list[str] markers: Genes.
    :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
    """
    vc =  VitessceConfig(schema_version="1.0.15", name='heatmap')
    adata = adata
    groupby = groupby

    if "markers" not in kwargs.keys():
        markers = []
    
    if type(kwargs.get("markers")) == str:
        markers = [kwargs.get("markers", [])]
    elif type(kwargs.get("markers")) == list: 
        markers = kwargs.get("markers", [])
        
    zarr_filepath = _create_zarr_filepath(adata, "violin")

    dataset = vc.add_dataset(name='data').add_object(AnnDataWrapper(
            adata_path=zarr_filepath,
            obs_set_paths=[f"obs/{groupby}"],
            obs_set_names=["cell type"],
            obs_feature_matrix_path="X"
        ))

    if type(markers) == list and len(markers) > 1:
        for gene in markers:
            genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset)
            cells = vc.add_view(cm.OBS_SETS, dataset=dataset)
            violin = vc.add_view('obsSetFeatureValueDistribution', dataset=dataset, uid=f'violin-plot-{gene}')
            vc.link_views(
            [violin, genes, cells], 
            ["featureSelection", "obsSetSelection"],
            [[gene], None]
            )
            vc.layout(hconcat(violin, genes, cells, split = [2,1,1]))
    else:
        genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset)
        cells = vc.add_view(cm.OBS_SETS, dataset=dataset)
        violin = vc.add_view('obsSetFeatureValueDistribution', dataset=dataset, uid='violin-plot')

        if "markers" in kwargs.keys():
            # print(markers)
            vc.link_views(
            [violin, genes, cells], 
            ["featureSelection"],
            [markers]
            )
    
        vc.layout(violin | genes / cells)

    vw = vc.widget()
    return vw

def dotplot(adata, groupby, **kwargs):
        """
        Creates interactive dotplot.

        :param AnnData adata: AnnData object.
        :param str groupby: Category group.
        :param list[str] markers: List of genes.
        :returns: Vitessce widget. Documentation can be found `here. <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ 
        """
        adata = adata
        groupby = groupby
    
        if "markers" in kwargs.keys():
            markers = kwargs["markers"]
    
        vc = VitessceConfig(schema_version="1.0.17", name='dotplot data')
    
        zarr_filepath = _create_zarr_filepath(adata, "dotplot")
    
        dataset = vc.add_dataset('dotplot data').add_object(AnnDataWrapper(
                adata_path = zarr_filepath,
                obs_set_paths=[f"obs/{groupby}"],
                obs_set_names=["cell type"],
                #obs_embedding_paths=["obsm/X_umap"],
                #obs_embedding_names=[""],
                obs_feature_matrix_path="X" 
    )).add_object(AnnDataWrapper(
        adata_path=zarr_filepath,
        coordination_values={
        "obsType": 'cell',
        "featureType": 'gene',
        "featureValueType": 'value',
        "sampleType": 'sample',
        }
        ))
    
    
        obsSets = vc.add_view(cm.OBS_SETS, dataset=dataset)
        
        featureList = vc.add_view(cm.FEATURE_LIST, dataset=dataset).set_props(enableMultiSelect=True)
        #violinPlots = vc.add_view('obsSetFeatureValueDistribution', dataset=dataset, uid='violin-plot')
        dotPlot = vc.add_view('dotPlot', dataset=dataset, uid='dot-plot')
        
        if "markers" in kwargs.keys():
            vc.link_views(
            [dotPlot, featureList], 
            ["featureSelection"],
            [markers]
        )
        
        
        vc.layout(dotPlot | featureList / (obsSets))
        vw = vc.widget()
        return vw

def configure_plots(disable_plots=[], enable_plots=[]): 
    """
    Deactivates and reactivates interactive Vitessce plots.

    :param list disable_plots: List of plots.
    :param list enable_plots: List of plots.
    """
    if any(plot in enable_plots for plot in disable_plots):
            raise RuntimeError("Plots cannot be in enable_plots and disable_plots simultaneously.")
        
    def monkeypatch(cls, func):
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

    def undo_monkeypatch(cls, func_name):
        orig_func_name = f"_orig_{func_name}"
        if hasattr(cls, orig_func_name):
            orig_func = getattr(cls, orig_func_name)
            setattr(cls, func_name, orig_func)

    def monkeypatch_spatialdata():
        orig_pl_name = "_orig_pl" # for storing pl
        if not hasattr(SpatialData, orig_pl_name):
            setattr(SpatialData, orig_pl_name, _CachedAccessor('pl', SpatialData.pl))
        setattr(SpatialData, 'pl', _CachedAccessor('pl', VitessceSpatialData))
    
    def undo_monkeypatch_spatialdata():
        orig_pl_name = "_orig_pl"
        if hasattr(SpatialData, '_orig_pl'):
            setattr(SpatialData, 'pl', _CachedAccessor('pl', SpatialData.pl))

    enable_embedding = True
    enable_spatial = True
    enable_dotplot = True
    enable_heatmap = True
    enable_violin = True
    
    enable_embedding = not "embedding" in disable_plots or "embedding" in enable_plots
    enable_spatial = not "spatial" in disable_plots or "spatial" in enable_plots
    enable_dotplot = not "dotplot" in disable_plots or "dotplot" in enable_plots
    enable_heatmap = not "heatmap" in disable_plots or "heatmap" in enable_plots
    enable_violin = not "violin" in disable_plots or "violin" in enable_plots
    
        
    if enable_embedding:
        monkeypatch(sc.pl, embedding)
    else:
        undo_monkeypatch(sc.pl, "embedding")
        print("deactivated Vitessce embedding")
        
    if enable_spatial:
        monkeypatch(sc.pl, spatial)
        monkeypatch_spatialdata()
    else:
        undo_monkeypatch(sc.pl, "spatial")
        undo_monkeypatch_spatialdata()
        print("deactivated Vitessce spatial")
        
    if enable_dotplot:
        monkeypatch(sc.pl, dotplot)
    else:
        undo_monkeypatch(sc.pl, "dotplot")
        print("deactivated Vitessce dotplot")
        
    if enable_heatmap:
        monkeypatch(sc.pl, heatmap)
    else:
        undo_monkeypatch(sc.pl, "heatmap")
        print("deactivated Vitessce heatmap")
        
    if enable_violin:
        monkeypatch(sc.pl, violin)
    else:
        undo_monkeypatch(sc.pl, "violin")
        print("deactivated Vitessce violin")
