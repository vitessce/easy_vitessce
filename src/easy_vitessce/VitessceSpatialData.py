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

from os.path import join

from spatialdata_plot.pl.basic import PlotAccessor

class VitessceSpatialData:
    """
    A class for configuring spatial plot with similar syntax to spatialdata from scverse.
    """

    # This is a class variable to determine whether the monkeypatching is enabled.
    # This is a workaround since our monkeypatching does not work with the existing instances of the SpatialData class.
    # In other words, when we change SpatialData.pl, the existing instances of SpatialData class are not affected.
    # Instead, we use this class variable.
    # This way, existing instances of the SpatialData class in which SpatialData.pl has been monkeypatched with VitessceSpatialData,
    # will see that monkeypatching is enabled/disabled, and will behave accordingly.
    _is_enabled = True

    def __init__(self, sdata):
        """
        Initializes filepaths, Vitessce configuration, SpatialDataWrapper arguments, color map.

        :param str spatialdata_filepath: filepath of spatialdata zarr file containing image data.
        :returns: Self, allows for chaining.
        """
        self.sdata = sdata
        if sdata.is_backed() and sdata.is_self_contained():
            self.sdata_filepath = sdata.path
        else:
            self.sdata_filepath = join("data", "sdata.zarr")
            sdata.write(self.sdata_filepath, overwrite=True)
        
        self.color = ""
        self.kwargs = {"sdata_path": self.sdata_filepath,
                # The following paths are relative to the root of the SpatialData zarr store on-disk.
                "table_path":"tables/table",
                "obs_feature_matrix_path":"tables/table/X",
                "coordinate_system":"global",
                "coordination_values":{
                    # The following tells Vitessce to consider each observation as a "spot"
                    "obsType": "spot",
                }}
        self.views = {"featureValueColormap": "viridis"} 

        # This is the static PlotAccessor instance that will be used when monkeypatching is not enabled.
        self._pl = PlotAccessor(sdata)
    
        
    def render_images(self, element="", **kwargs):
        """
        Renders image.

        :param str element: location of image data inside "images" folder.
        :returns: Self, allows for chaining.
        """
        if not VitessceSpatialData._is_enabled:
            return self._pl.render_images(element=element, **kwargs)

        self.image = f"images/{element}"
        self.image_path = {"image_path":f"images/{element}"}
        self.kwargs.update(self.image_path)

        return self.sdata
        
    def render_shapes(self, element="", **kwargs):
        """
        Renders shapes, e.g. "cells".

        :param str element: location of shape data inside "shapes" folder.
        :param str color: gene.
        :param str cmap: color map (viridis, plasma, jet).
        :returns: Self, allows for chaining.
        """
        if not VitessceSpatialData._is_enabled:
            return self._pl.render_shapes(element=element, **kwargs)

        self.color = kwargs.get("color", "")

        if (self.sdata.shapes[element]["geometry"].geom_type.iloc[0]) == 'Polygon': # vitessce only has polygon and circles
            obs_path = {"obs_segmentations_path": f"shapes/{element}"}
            segmentation_layer = {"spatial_segmentation_layer": [obs_path, self.image]}
            self.kwargs.update(segmentation_layer)
        else:
            obs_path = {"obs_spots_path": f"shapes/{element}"}
            
        self.kwargs.update(obs_path)
        print(f"self.kwargs: {self.kwargs}")

        if "table_layer" in kwargs.keys():
            # have user specify which file to use?
            table_path = {"table_path": f"tables/{kwargs.get('table_layer')}"}
            matrix_path = {"obs_feature_matrix_path": f"tables/{kwargs.get('table_layer')}/X"}
            self.kwargs.update(table_path)
            self.kwargs.update(matrix_path)

        if "color" in kwargs.keys():
            self.color = kwargs.get("color")
            if self.color in self.sdata.tables["table"].var.index: # gene
                color = {"featureSelection": [kwargs["color"]]}
                color_encoding = {"obsColorEncoding": "geneSelection"}
                
                self.views.update(color)
                self.views.update(color_encoding)
                print(f"self.views: {self.views}")
                
            elif self.color in self.sdata.tables["table"].obs: # categorical?
                color = {"obsSetSelection": [[kwargs["color"]]]}
                color_encoding = {"obsColorEncoding": "cellSetSelection"}
                
                self.views.update(color)
                self.views.update(color_encoding)
                print(f"self.views: {self.views}")
            
        if "cmap" in kwargs.keys():
            cmap = {"featureValueColormap": kwargs["cmap"]}
            self.views.update(cmap)
            
        return self.sdata

    def render_labels(self, element="", **kwargs):
        """
        Renders label data.

        :param str element: location of label data in "labels" folder.
        :returns: Self, allows for chaining.
        """
        if not VitessceSpatialData._is_enabled:
            return self._pl.render_labels(element=element, **kwargs)
        
        labels_path = {"obs_segmentations_path":f"labels/{element}"}
        self.kwargs.update(labels_path)

        return self.sdata

    def render_points(self, element="", **kwargs):
        """
        Renders points.

        :param str element: location of point data in "points" folder.
        :returns: Self, allows for chaining.
        """
        if not VitessceSpatialData._is_enabled:
            return self._pl.render_points(element=element, **kwargs)
        
        obs_points_path = {"obs_points_path":f"points/{element}"}
        self.kwargs.update(obs_points_path)

        return self.sdata
    
    def show(self, **kwargs):
        """
        Displays spatial plot.
        
        :returns: Vitessce widget.
        """
        if not VitessceSpatialData._is_enabled:
            return self._pl.show(**kwargs)
            
        self.vc = VitessceConfig(schema_version="1.0.18", name='spatial data')
        self.wrapper = SpatialDataWrapper(**self.kwargs)
        
        dataset = self.vc.add_dataset(name='Spatial Data').add_object(self.wrapper)
        
        side_list = vt.OBS_SETS if self.color in self.sdata.tables["table"].obs else vt.FEATURE_LIST
        print(self.color)

        # Add views (visualizations) to the configuration:
        spatial = self.vc.add_view("spatialBeta", dataset=dataset)
        feature_obs_list = self.vc.add_view(side_list, dataset=dataset)
        layer_controller = self.vc.add_view("layerControllerBeta", dataset=dataset)
        
        self.vc.link_views([spatial, layer_controller, feature_obs_list], ['obsType'], [self.wrapper.obs_type_label])
        
        self.vc.link_views_by_dict([spatial, layer_controller, feature_obs_list], self.views, meta=False)
        
        # Layout the views
        self.vc.layout(spatial | (feature_obs_list / layer_controller))
        
        return self.vc.widget()