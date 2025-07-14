from vitessce import (
    VitessceConfig,
    AnnDataWrapper,
    SpatialDataWrapper,
    ImageOmeZarrWrapper,
    CoordinationLevel as CL,
    Component as cm,
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

from anndata import (AnnData, read_h5ad)


def configure_plots(enable_plots=[], disable_plots=[]):
    
    def monkeypatch(cls):
        def decorator(func):
            func_name = func.__name__
            orig_func_name = f"_orig_{func_name}"
            if not hasattr(cls, orig_func_name):
              orig_func = getattr(cls, func_name)
              setattr(cls, orig_func_name, orig_func)
            setattr(cls, func_name, func)
            return func
        return decorator

    def undo_monkeypatch(cls, func_name):
        orig_func_name = f"_orig_{func_name}"
        if hasattr(cls, orig_func_name):
            orig_func = getattr(cls, orig_func_name)
           # print(f"func_name: {func_name}")
            #print(f"orig_func: {orig_func}")
            setattr(cls, func_name, orig_func)
            #print("hasattr check")
            
    enable_embedding = True
    enable_spatial = True
    enable_dotplot = True
    enable_heatmap = True
    enable_violin = True

    if "embedding" in disable_plots:
        enable_embedding = False
    elif "embedding" in enable_plots:
        enable_embedding = True

    if "spatial" in disable_plots:
        enable_spatial = False
    elif "spatial" in enable_plots:
        enable_spatial = True

    if "dotplot" in disable_plots:
        enable_dotplot = False
    elif "dotplot" in enable_plots:
        enable_dotplot = True
        
    if "heatmap" in disable_plots:
        enable_heatmap = False
    elif "heatmap" in enable_plots:
        enable_heatmap = True

    if "violinl" in disable_plots:
        enable_violin = False
    elif "violin" in enable_plots:
        enable_violin = True

            
    if not enable_embedding:
        undo_monkeypatch(sc.pl, "embedding")
        print("disabled Vitessce embedding")

                   
    if enable_embedding:
        # print("doesn't run?")
        @monkeypatch(sc.pl)
        def embedding(adata, basis, **kwargs):
            basis = basis
            adata = adata
        
            color = kwargs.get("color", "")
            color_map = kwargs.get("color_map", "viridis")
            size  = kwargs.get("size", 2.5)
        
            if "color_map" in kwargs.keys():
                if "plasma" not in kwargs["color_map"].lower() and "viridis" not in kwargs["color_map"].lower() and "jet" not in kwargs["color_map"].lower():
                    print("Invalid color_map. Supported color_maps: plasma, viridis, jet. Default set to plasma.")
                    color_map = "plasma"
            
            def create_zarr_filepath(adata): 
                    zarr_filepath = join("data", "embedding_file.zarr")
                    if os.path.exists(zarr_filepath) and os.path.isdir(zarr_filepath):
                        shutil.rmtree(zarr_filepath)   
                    adata.write_zarr(zarr_filepath, chunks=[adata.shape[0], VAR_CHUNK_SIZE])
                    return zarr_filepath
        
            zarr_filepath = create_zarr_filepath(adata)
        
            if basis == "umap" or basis == "pca":
                vc =  VitessceConfig(schema_version="1.0.15", name='UMAP' if basis=="umap" else "PCA")
                
                dataset = vc.add_dataset(name='data').add_object(AnnDataWrapper(
                    adata_store=zarr_filepath,
                    obs_embedding_paths=["obsm/X_umap" if basis == "umap" else "obsm/X_pca"],
                    obs_embedding_names=["UMAP" if basis == "umap" else "PCA"],
                    obs_feature_matrix_path="X"
                ))
           
                mapping = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="UMAP" if basis=="umap" else "PCA") # mapping value corresponds to one of the obs_embedding_names values.
                genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset)
                    
                if ("heatmap" in kwargs) and (kwargs["heatmap"] == True):
                    heatmap = vc.add_view(cm.HEATMAP, dataset=dataset)
                    vc.link_views(
                        [mapping, genes, heatmap], 
                        ["featureSelection", "obsColorEncoding", "embeddingObsRadiusMode", "embeddingObsRadius", "featureValueColormap"], # https://vitessce.io/docs/coordination-types/
                        [[color], "geneSelection", "manual", size , color_map]
                    )
                    vc.layout(mapping | genes / heatmap);
                
                else:
                    vc.link_views(
                        [mapping, genes], 
                        ["featureSelection", "obsColorEncoding", "embeddingObsRadiusMode", "embeddingObsRadius", "featureValueColormap"], # https://vitessce.io/docs/coordination-types/
                        [[color], "geneSelection", "manual", size , color_map]
                    )
                    vc.layout(mapping | genes);
                   
                vw = vc.widget()
                return vw
        
        
            elif basis == "tsne":
                if "X_tsne" not in adata.obsm:
                    adata = sc.tl.tsne(adata, random_state=1) # add a check for this in case it wasn't ran previously
               
                vc =  VitessceConfig(schema_version="1.0.15", name='t-SNE')
                
                if "color" in kwargs and kwargs["color"] in adata.obs.columns:
                # if "color" in kwargs and kwargs["color"] in adata.obs.columns:
                    dataset = vc.add_dataset(name='tsne data').add_object(AnnDataWrapper(
                            adata_store=zarr_filepath,
                            obs_set_paths=[f"obs/{kwargs['color']}"],
                            obs_set_names=["colors"],
                            obs_embedding_paths=["obsm/X_tsne"],
                            obs_embedding_names=["t-SNE"],
                            obs_feature_matrix_path="X"
                        ))
                    cells = vc.add_view(cm.OBS_SETS, dataset=dataset)
                
                else:
                        dataset = vc.add_dataset(name='tsne data').add_object(AnnDataWrapper(
                        adata_store=zarr_filepath,
                        # obs_set_paths=["obs/bulk_labels"], # consider the order of kwargs
                        # obs_set_names=["colors"],
                        obs_embedding_paths=["obsm/X_tsne"],
                        obs_embedding_names=["t-SNE"],
                        obs_feature_matrix_path="X"
                    ))
                        cells = vc.add_view(cm.OBS_SETS, dataset=dataset)
                        if len(color) != 0:
                            print("Color label not found in data")
               
                tsne = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="t-SNE") # mapping value corresponds to one of the obs_embedding_names values.
                # genes = vc.add_view(cm.OBS_SETS, dataset=dataset)
                
                # https://python-docs.vitessce.io/api_config.html?highlight=link_views#vitessce.config.VitessceConfig.link_views 
                vc.link_views(
                    [tsne, cells], 
                    ["obsColorEncoding", "embeddingObsRadiusMode", "embeddingObsRadius", "embeddingObsOpacityMode", "embeddingObsOpacity"], # https://vitessce.io/docs/coordination-types/
                    ["cellSetSelection", "manual", size, "manual", 0.6]
                )
                
                vc.layout(tsne | cells);
                
                vw = vc.widget(port=9000)
                return vw
            
    if not enable_spatial:
        undo_monkeypatch(sc.pl, "spatial")
        print("disabled Vitessce spatial")
    
    elif enable_spatial:
            @monkeypatch(sc.pl)
            def spatial(adata, **kwargs):
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
        
    if not enable_dotplot:
        undo_monkeypatch(sc.pl, "dotplot")
        print("disabled Vitessce dotplot")
    
    elif enable_dotplot:
            @monkeypatch(sc.pl)
            def dotplot(adata, **kwargs):
                adata = adata
            
                if "markers" in kwargs.keys():
                    markers = kwargs["markers"]
                
                if "groupby" in kwargs.keys():
                    groupby = kwargs["groupby"]
            
                vc = VitessceConfig(schema_version="1.0.17", name='dotplot data')
            
                def create_zarr_filepath(adata): 
                        zarr_filepath = join("data", "dotplot_file.zarr")
                        if os.path.exists(zarr_filepath) and os.path.isdir(zarr_filepath):
                            shutil.rmtree(zarr_filepath)   
                        adata.write_zarr(zarr_filepath, chunks=[adata.shape[0], VAR_CHUNK_SIZE])
                        return zarr_filepath
            
                zarr_filepath = create_zarr_filepath(adata)
            
                dataset = vc.add_dataset('dotplot data').add_object(AnnDataWrapper(
                        adata_path = zarr_filepath,
                        obs_set_paths=[f"obs/{groupby}"],
                        obs_set_names=["cell type"],
                        #obs_embedding_paths=["obsm/X_umap"],
                        #obs_embedding_names=[""],
                        obs_feature_matrix_path="X" # "layers/logcounts" doesn't exist
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
                
                if "markers" in kwargs:
                    vc.link_views(
                    [dotPlot, featureList], 
                    ["featureSelection"],
                    [markers]
                )
                
                
                vc.layout(dotPlot | featureList / (obsSets))
                vw = vc.widget()
                return vw

    if not enable_heatmap:
        undo_monkeypatch(sc.pl, "heatmap")
        print("disabled Vitessce heatmap")
    
    elif enable_heatmap:
        @monkeypatch(sc.pl)
        def heatmap(adata, **kwargs):
            vc =  VitessceConfig(schema_version="1.0.15", name='heatmap')
            adata = adata
            color_map = kwargs.get("color_map", "viridis")

            if "groupby" in kwargs:
                groupby = kwargs["groupby"]

            if "markers" in kwargs:
                markers = kwargs["markers"]
                adata.var["genes"] = list(adata.var.index)
                adata.var["in_markers"] = adata.var["genes"].apply(lambda gene: True if gene in markers else False)
                
            def create_zarr_filepath(adata): 
                    zarr_filepath = join("data", "heatmap_file.zarr")
                    if os.path.exists(zarr_filepath) and os.path.isdir(zarr_filepath):
                        shutil.rmtree(zarr_filepath)   
                    adata.write_zarr(zarr_filepath, chunks=[adata.shape[0], VAR_CHUNK_SIZE])
                    return zarr_filepath

            zarr_filepath = create_zarr_filepath(adata)

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
                    
            vc.layout(heatmap | cells);

            vw = vc.widget()
            return vw
        
    if not enable_violin:
         undo_monkeypatch(sc.pl, "violin")
         print("disabled Vitessce violin plot")

    elif enable_violin:
        @monkeypatch(sc.pl)
        def violin(adata, **kwargs):
            vc =  VitessceConfig(schema_version="1.0.15", name='heatmap')
            adata = adata

            if "markers" in kwargs:
                markers = kwargs["markers"]
                
            def create_zarr_filepath(adata): 
                    zarr_filepath = join("data", "violin_file.zarr")
                    if os.path.exists(zarr_filepath) and os.path.isdir(zarr_filepath):
                        shutil.rmtree(zarr_filepath)   
                    adata.write_zarr(zarr_filepath, chunks=[adata.shape[0], VAR_CHUNK_SIZE])
                    return zarr_filepath

            zarr_filepath = create_zarr_filepath(adata)

            dataset = vc.add_dataset(name='data').add_object(AnnDataWrapper(
                    adata_path=zarr_filepath,
                    obs_set_paths=["obs/bulk_labels"],
                    obs_set_names=["cell type"],
                    obs_feature_matrix_path="X"
                ))

            genes = vc.add_view(cm.FEATURE_LIST, dataset=dataset).set_props(enableMultiSelect=True)
            cells = vc.add_view(cm.OBS_SETS, dataset=dataset)
            violin = vc.add_view('obsSetFeatureValueDistribution', dataset=dataset, uid='violin-plot')

            if "markers" in kwargs:
                vc.link_views(
                [violin, genes, cells], 
                ["featureSelection"],
                [markers]
                )
            
            vc.layout(violin | genes / cells);

            vw = vc.widget()
            return vw
        
                            