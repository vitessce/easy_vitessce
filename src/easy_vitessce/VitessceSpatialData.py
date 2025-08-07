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

from os import join

class VitessceSpatialData:
    """
    A class for configuring spatial plot with similar syntax to spatialdata from scverse.
    """
    def __init__(self, sdata):
        """
        Initializes filepaths, Vitessce configuration, SpatialDataWrapper arguments, color map.

        :param str spatialdata_filepath: filepath of spatialdata zarr file containing image data.
        :returns: Self, allows for chaining.
        """
        self.sdata_filepath = join("data", "sdata.zarr")
        sdata.write(self.sdata_filepath, overwrite=True)
        self.vc = VitessceConfig(schema_version="1.0.18", name='spatial data')
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
        """
        Renders image.

        :param str element: location of image data inside "images" folder.
        :returns: Self, allows for chaining.
        """
        image_path = {"image_path":f"images/{element}"}
        self.kwargs.update(image_path)

        return self
        
    def render_shapes(self, element="", **kwargs):
        """
        Renders shapes, e.g. "cells".

        :param str element: location of shape data inside "shapes" folder.
        :param str color: gene.
        :param str color_map: color map (viridis, plasma, jet).
        :returns: Self, allows for chaining.
        """
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
        """
        Renders label data.

        :param str element: location of label data in "labels" folder.
        :returns: Self, allows for chaining.
        """
        labels_path = {"obs_segmentations_path":f"labels/{element}"}
        self.kwargs.update(labels_path)

        return self
    
    def show(self):
        """
        Displays spatial plot.
        
        :returns: Vitessce widget.
        """
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