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
    CoordinationLevel as CL,
    SpatialDataWrapper,
    get_initial_coordination_scope_prefix,
)

from os.path import join

from spatialdata_plot.pl.basic import PlotAccessor
from spatialdata import get_element_annotators

from easy_vitessce.widget import _to_widget, config

# This class is analogous to PlotAccessor from spatialdata-plot.
# Reference: https://github.com/scverse/spatialdata-plot/blob/788eb2206cca8f4c21977c4f7b08a818ee6580f7/src/spatialdata_plot/pl/basic.py#L68
class VitesscePlotAccessor:
    """
    A class for configuring a spatial plot, using the same syntax as spatialdata-plot.
    """

    # This is a class variable to determine whether the monkeypatching is enabled.
    # This is a workaround since our monkeypatching does not work with the existing instances of the SpatialData class.
    # In other words, when we change SpatialData.pl, the existing instances of SpatialData class are not affected.
    # Instead, we use this class variable.
    # This way, existing instances of the SpatialData class in which SpatialData.pl has been monkeypatched with VitesscePlotAccessor,
    # will see that monkeypatching is enabled/disabled, and will behave accordingly.
    _is_enabled = True

    def __init__(self, sdata):
        """
        Initialize the plot accessor.

        :param SpatialData sdata: The SpatialData object to use for plotting.
        """
        self.sdata = sdata
        if sdata.is_backed() and sdata.is_self_contained():
            self.sdata_filepath = sdata.path
        else:
            self.sdata_filepath = join(config.get('data.out_dir'), "sdata.zarr")
            sdata.write(self.sdata_filepath, overwrite=config.get('data.overwrite'))
        
        self._init_params()

        # This is the static PlotAccessor instance that will be used when monkeypatching is not enabled.
        self._pl = PlotAccessor(sdata)
    
    def _init_params(self):
        self.obs_type = "cell" # TODO: support multiple obs types (one per layer?)
        self.wrapper_args = {
            "sdata_path": self.sdata_filepath,
            # The following paths are relative to the root of the SpatialData zarr store on-disk.
            "table_path":"tables/table",
            "obs_feature_matrix_path":"tables/table/X",
            "coordinate_system":"global",
            "coordination_values":{}
        }

        self.has_gene_color_encoding = False
        self.has_cellset_color_encoding = False

        self.global_coordination = {"featureValueColormap": "viridis", "obsColorEncoding": "geneSelection"}

        self.image_layer_coordination = []
        self.segmentation_layer_coordination = []
        self.spot_layer_coordination = []
        self.point_layer_coordination = []
        
    def render_images(self, element=None, **kwargs):
        """
        Renders image.

        :param str element: location of image data inside "images" folder.
        :returns: Self, allows for chaining.
        """
        if not VitesscePlotAccessor._is_enabled:
            return self._pl.render_images(element=element, **kwargs)

        # channel (list[str] | list[int] | str | int | None)
        #   To select specific channels to plot.
        #   Can be a single channel name/int or a list of channel names/ints.
        #   If None, all channels will be used.
        channel_param = kwargs.get("channel", None)
        # cmap (list[Colormap | str] | Colormap | str | None)
        #   Colormap or list of colormaps for continuous annotations, see matplotlib.colors.Colormap.
        #   Each colormap applies to a corresponding channel.
        cmap_param = kwargs.get("cmap", None)
        # palette (list[str] | str | None)
        #   Palette to color images.
        #   The number of palettes should be equal to the number of channels.
        palette_param = kwargs.get("palette", None)
        # alpha (float | int, default 1.0)
        #   Alpha value for the images.
        #   Must be a numeric between 0 and 1.
        alpha_param = kwargs.get("alpha", None)

        self.image = f"images/{element}"
        self.image_path = {"image_path":f"images/{element}"}
        self.wrapper_args.update(self.image_path)

        # Palette logic in spatialdata-plot:
        # Reference: https://github.com/scverse/spatialdata-plot/blob/010560f7eebdd245693a8c55eede0f895a636f5c/src/spatialdata_plot/pl/utils.py#L685

        # RGB vs. non-RGB logic in spatialdata-plot:
        # Reference: https://github.com/scverse/spatialdata-plot/blob/010560f7eebdd245693a8c55eede0f895a636f5c/src/spatialdata_plot/pl/render.py#L865
        img = self.sdata.images[element]
        channels = img.coords["c"].values.tolist() if channel_param is None else channel_param

        # the channel parameter has been previously validated, so when not None, render_params.channel is a list
        assert isinstance(channels, list)
        n_channels = len(channels)

        # Not ideal logic. Should ideally only use the OME-NGFF color model metadata. But this is what spatialdata-plot does.
        photometric_interpretation = "RGB" if palette_param is None and channel_param is None and n_channels == 3 else "BlackIsZero"

        self.image_layer_coordination = [
            # We want to keep any existing spatial layer coordination information.
            *self.image_layer_coordination,
            {
                'spatialLayerOpacity': alpha_param if alpha_param is not None else 1.0,
                'photometricInterpretation': photometric_interpretation,
                # 'imageChannel': [{
                #     # TODO: specify spatialTargetC if channel_param is not None
                #     'spatialChannelColor': [255, 255, 255], # TODO: use the palette or cmap
                # }]
            },
        ]

        return self.sdata
        
    def render_shapes(self, element="", **kwargs):
        """
        Renders shapes, e.g. "cells".

        :param str element: location of shape data inside "shapes" folder.
        :param str color: gene.
        :param str cmap: color map (viridis, plasma, jet).
        :returns: Self, allows for chaining.
        """
        if not VitesscePlotAccessor._is_enabled:
            return self._pl.render_shapes(element=element, **kwargs)
        
        color_param = kwargs.get("color")

        # vitessce only has polygon and circles
        if self.sdata.shapes[element]["geometry"].geom_type.iloc[0] == 'Polygon':
            # This is a polygon-type Shapes element, so we use obs_segmentations_path.
            obs_path = {"obs_segmentations_path": f"shapes/{element}"}

            self.segmentation_layer_coordination = [
                # We want to keep any existing spatial layer coordination information.
                *self.segmentation_layer_coordination,
                {
                    'segmentationChannel': [{
                        # We initialize with a single channel.

                    }],
                },
            ]
        else:
            self.obs_type = "spot"
            # This is a circle-type Shapes element, so we use obs_spots_path.
            obs_path = {"obs_spots_path": f"shapes/{element}"}

            self.spot_layer_coordination = [
                # We want to keep any existing spatial layer coordination information.
                *self.spot_layer_coordination,
                {
                    
                },
            ]
        
        self.wrapper_args.update(obs_path)

        table_name = kwargs.get("table_name", None)
        if table_name is None:
            annotating_tables = list(get_element_annotators(self.sdata, element))
            if len(annotating_tables) > 0:
                # Use the first annotating table if no specific table is provided.
                table_name = annotating_tables[0]

        if table_name is not None:
            # have user specify which matrix to use?
            table_path = {"table_path": f"tables/{table_name}"}
            self.wrapper_args.update(table_path)

            self.wrapper_args = {
                **self.wrapper_args,
                # TODO: check for X existence first?
                "obs_feature_matrix_path": f"tables/{table_name}/X"
            }

            # TODO: configure all obsSets in the table here, to allow the user to select them regardless of the "color" parameter value,
            # rather than only when the "color" parameter is set to a categorical obs column (down below).

        if color_param is not None:
            if table_name is None:
                raise ValueError("The 'color' parameter was provided, but an annotating table was not found. You may need to specify 'table_name' explicitly.")
            
            if color_param in self.sdata.tables[table_name].var.index: # gene
                self.has_gene_color_encoding = True
                color = {"featureSelection": [color_param]}
                color_encoding = {"obsColorEncoding": "geneSelection"}
                
                self.global_coordination.update(color)
                self.global_coordination.update(color_encoding)

            elif color_param in self.sdata.tables[table_name].obs: # categorical?
                self.has_cellset_color_encoding = True
                # TODO: depends on https://github.com/vitessce/vitessce/issues/2254
                color = {"obsSetSelection": [[color_param]]}
                color_encoding = {"obsColorEncoding": "cellSetSelection"}
                
                self.global_coordination.update(color)
                self.global_coordination.update(color_encoding)

                # Here we configure obsSets for self.wrapper_args
                self.wrapper_args = {
                    **self.wrapper_args,
                    "obs_set_paths": [f"tables/{table_name}/obs/{color_param}"],
                    "obs_set_names": [color_param],
                }
            else:
                # TODO: support a static color, such as "red" or "#FF0000"?
                raise ValueError(f"Color value did not map to a value in var.index or obs.columns of table {table_name}.")
            
        if "cmap" in kwargs.keys():
            cmap = {"featureValueColormap": kwargs["cmap"]}
            self.global_coordination.update(cmap)
            
        return self.sdata

    def render_labels(self, element="", **kwargs):
        """
        Renders label data.

        :param str element: location of label data in "labels" folder.
        :returns: Self, allows for chaining.
        """
        if not VitesscePlotAccessor._is_enabled:
            return self._pl.render_labels(element=element, **kwargs)
        
        labels_path = {"obs_segmentations_path":f"labels/{element}"}
        self.wrapper_args.update(labels_path)

        # Same coloring logic as in render_shapes.
        color_param = kwargs.get("color")

        self.segmentation_layer_coordination = [
            # We want to keep any existing spatial layer coordination information.
            *self.segmentation_layer_coordination,
            {
                'segmentationChannel': [{
                    # We initialize with a single channel.

                }],
            },
        ]

        table_name = kwargs.get("table_name", None)
        if table_name is None:
            annotating_tables = list(get_element_annotators(self.sdata, element))
            if len(annotating_tables) > 0:
                # Use the first annotating table if no specific table is provided.
                table_name = annotating_tables[0]

        if table_name is not None:
            # have user specify which matrix to use?
            table_path = {"table_path": f"tables/{table_name}"}
            self.wrapper_args.update(table_path)

            self.wrapper_args = {
                **self.wrapper_args,
                # TODO: check for X existence first?
                "obs_feature_matrix_path": f"tables/{table_name}/X"
            }

            # TODO: configure all obsSets in the table here, to allow the user to select them regardless of the "color" parameter value,
            # rather than only when the "color" parameter is set to a categorical obs column (down below).

        if color_param is not None:
            if table_name is None:
                raise ValueError("The 'color' parameter was provided, but an annotating table was not found. You may need to specify 'table_name' explicitly.")
            
            if color_param in self.sdata.tables[table_name].var.index: # gene
                self.has_gene_color_encoding = True
                color = {"featureSelection": [color_param]}
                color_encoding = {"obsColorEncoding": "geneSelection"}
                
                self.global_coordination.update(color)
                self.global_coordination.update(color_encoding)

            elif color_param in self.sdata.tables[table_name].obs: # categorical?
                self.has_cellset_color_encoding = True
                # TODO: depends on https://github.com/vitessce/vitessce/issues/2254
                color = {"obsSetSelection": [[color_param]]}
                color_encoding = {"obsColorEncoding": "cellSetSelection"}
                
                self.global_coordination.update(color)
                self.global_coordination.update(color_encoding)

                # Here we configure obsSets for self.wrapper_args
                self.wrapper_args = {
                    **self.wrapper_args,
                    "obs_set_paths": [f"tables/{table_name}/obs/{color_param}"],
                    "obs_set_names": [color_param],
                }
            else:
                # TODO: support a static color, such as "red" or "#FF0000"?
                raise ValueError(f"Color value did not map to a value in var.index or obs.columns of table {table_name}.")
            
        if "cmap" in kwargs.keys():
            cmap = {"featureValueColormap": kwargs["cmap"]}
            self.global_coordination.update(cmap)
            
        return self.sdata

    def render_points(self, element="", **kwargs):
        """
        Renders points.

        :param str element: location of point data in "points" folder.
        :returns: Self, allows for chaining.
        """
        if not VitesscePlotAccessor._is_enabled:
            return self._pl.render_points(element=element, **kwargs)
        
        obs_points_path = {"obs_points_path":f"points/{element}"}
        self.wrapper_args.update(obs_points_path)

        self.point_layer_coordination = [
            # We want to keep any existing spatial layer coordination information.
            *self.point_layer_coordination,
            {
                
            },
        ]

        return self.sdata
    
    def show(self, coordinate_systems=None, **kwargs):
        """
        Displays spatial plot.
        
        :returns: Vitessce widget. Learn more at the vitessce-python `docs <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ .
        """
        if not VitesscePlotAccessor._is_enabled:
            return self._pl.show(**kwargs)
            
        self.vc = VitessceConfig(schema_version="1.0.18", name='spatial data')

        if not (coordinate_systems is None or isinstance(coordinate_systems, str)):
            raise NotImplementedError("A list of multiple 'coordinate_systems' is not yet supported.")

        self.wrapper = SpatialDataWrapper(**{
            **self.wrapper_args,
            **({ "coordinate_system": coordinate_systems } if coordinate_systems is not None else {}),
            "coordination_values": {
                **self.wrapper_args.get("coordination_values", {}),
                "obsType": self.obs_type,
            },
        })
        
        dataset_uid = "A"
        dataset = self.vc.add_dataset(name='Spatial Data', uid=dataset_uid).add_object(self.wrapper)

        side_list = vt.OBS_SETS if self.has_cellset_color_encoding else vt.FEATURE_LIST

        # Add views (visualizations) to the configuration:
        spatial = self.vc.add_view("spatialBeta", dataset=dataset)
        layer_controller = self.vc.add_view("layerControllerBeta", dataset=dataset)
        feature_obs_list = self.vc.add_view(side_list, dataset=dataset) # either feature list or obs sets # TODO: both obsSets and featureList

        spatial_views = [spatial, layer_controller]
        all_views = [spatial, layer_controller, feature_obs_list]

        # Create coordination scope objects for self.global_coordination
        ct_names = []
        ct_vals = []
        for ct_name, ct_val in self.global_coordination.items():
            ct_names.append(ct_name)
            ct_vals.append(ct_val)
        
        ct_scopes = self.vc.add_coordination(*ct_names)
        for i, ct_scope in enumerate(ct_scopes):
            ct_scope.set_value(ct_vals[i])

        # Link the views together
        self.vc.link_views(all_views, ['obsType'], [self.obs_type])
        self.vc.link_views_by_dict(all_views, dict(zip(ct_names, ct_scopes)), meta=True)
        self.vc.link_views_by_dict(spatial_views, {
            "imageLayer": CL(self.image_layer_coordination),
        }, meta=True, scope_prefix=get_initial_coordination_scope_prefix(dataset_uid, "image"))
        self.vc.link_views_by_dict(spatial_views, {
            "segmentationLayer": CL([
                {
                    **layer_dict,
                    'segmentationChannel': CL([
                        {
                            **channel_dict,
                            # TODO: limit this to the coordination types that are applicable to segmentationChannel objects
                            **dict(zip(ct_names, ct_scopes))
                        }
                        for channel_dict in layer_dict.get('segmentationChannel', [{}])
                    ])
                }
                for layer_dict in self.segmentation_layer_coordination
            ]),
        }, meta=True, scope_prefix=get_initial_coordination_scope_prefix(dataset_uid, "obsSegmentations"))
        self.vc.link_views_by_dict(spatial_views, {
            "spotLayer": CL(self.spot_layer_coordination),
        }, meta=True, scope_prefix=get_initial_coordination_scope_prefix(dataset_uid, "obsSpots"))
        self.vc.link_views_by_dict(spatial_views, {
            "pointLayer": CL(self.point_layer_coordination),
        }, meta=True, scope_prefix=get_initial_coordination_scope_prefix(dataset_uid, "obsPoints"))
        
        # Layout the views
        self.vc.layout(spatial | (feature_obs_list / layer_controller))
        
        vw = _to_widget(self.vc)

        # Cleanup
        self._init_params()

        return vw
