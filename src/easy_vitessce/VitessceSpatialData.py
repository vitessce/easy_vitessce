#import os
#from os.path import join, isfile, isdir
#from urllib.request import urlretrieve
#import zipfile
# import scanpy as sc
# import spatialdata as sd
# #import spatialdata_plot
# import numpy as np
# import matplotlib.pyplot as plt
# import shutil

from vitessce import (
    VitessceConfig,
    ViewType as vt,
    #CoordinationType as ct,
    #CoordinationLevel as CL,
    SpatialDataWrapper,
    #get_initial_coordination_scope_prefix
)

class VitessceSpatialData:
    def __init__(self, spatialdata_filepath, zip_filepath):
        self.sdata_filepath = spatialdata_filepath
        self.zip_filepath = zip_filepath
        self.vc = VitessceConfig(schema_version="1.0.16", name='Mouse Brain')
        self.kwargs = {"sdata_path": self.sdata_filepath,
                # The following paths are relative to the root of the SpatialData zarr store on-disk.
                "table_path":"tables/table",
                "obs_feature_matrix_path":"tables/table/X",
                "coordinate_syste":"global",
                "coordination_values":{
                    # The following tells Vitessce to consider each observation as a "spot"
                    "obsType": "spot",
                }}
        self.views = {"featureValueColormap": "viridis"}
        
        self.pl = self
        
        
    def render_images(self, element=""):
        image_path = {"image_path":f"images/{element}"}
        self.kwargs.update(image_path)

        return self
        
    def render_shapes(self, element="", **kwargs):
        obs_spots_path = {"obs_spots_path": f"shapes/{element}"}
        self.kwargs.update(obs_spots_path)

        if "color" in kwargs.keys():
            color = {"featureSelection": [kwargs["color"]]}
            color_encoding = {"obsColorEncoding": "geneSelection"}
            self.views.update(color)
            self.views.update(color_encoding)
            
        if "color_map" in kwargs.keys():
            color_map = {"featureValueColormap": kwargs["color_map"]}
            self.views.update(color_map)
            
        return self

    def render_labels(self, element=""):
        labels_path = {"labels_path":f"labels/{element}"}
        self.kwargs.update(labels_path)

        return self
    
    def show(self):
        self.wrapper = SpatialDataWrapper(**self.kwargs)
        
        dataset = self.vc.add_dataset(name='Mouse Brain Merfish').add_object(self.wrapper)
        
        # Add views (visualizations) to the configuration:
        spatial = self.vc.add_view("spatialBeta", dataset=dataset)
        feature_list = self.vc.add_view(vt.FEATURE_LIST, dataset=dataset)
        layer_controller = self.vc.add_view("layerControllerBeta", dataset=dataset)
        
        self.vc.link_views([spatial, layer_controller, feature_list], ['obsType'], [self.wrapper.obs_type_label])
        
        self.vc.link_views_by_dict([spatial, layer_controller, feature_list], self.views, meta=False)
        
        # Layout the views
        self.vc.layout(spatial | (feature_list / layer_controller))
        
        return self.vc.widget()