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
    hconcat,
    vconcat,
)

from os.path import join

from spatialdata_plot.pl.basic import PlotAccessor
from spatialdata import get_element_annotators

from easy_vitessce.widget import _to_widget, config
from easy_vitessce.colors import to_uint8_rgb

def shared_render_shapes_and_labels(
        sdata, element, table_name, table_layer, color, cmap, norm, groups, palette, obs_type, feature_type,
        # Note: These dict params are modified by this function.
        wrapper_args, obs_type_to_num_rows, feature_type_to_num_rows,
    ):

    if table_name is None:
        annotating_tables = list(get_element_annotators(sdata, element))
        if len(annotating_tables) > 0:
            # Use the first annotating table if no specific table is provided.
            table_name = annotating_tables[0]

    if table_name is not None:
        # have user specify which matrix to use?
        wrapper_args["table_path"] = f"tables/{table_name}"
        wrapper_args["obs_feature_matrix_path"] = f"tables/{table_name}/X" if table_layer is None else f"tables/{table_name}/layers/{table_layer}"

        obs_num_rows = sdata.tables[table_name].obs.shape[0]
        var_num_rows = sdata.tables[table_name].var.shape[0]

        if obs_type not in obs_type_to_num_rows:
            obs_type_to_num_rows[obs_type] = obs_num_rows
        elif obs_type_to_num_rows[obs_type] != obs_num_rows:
            # TODO: support automatically by using the element name as obsType?
            # Maybe force the user to configure something like ev.config.set({ "sdata_pl_use_element_name_for_entity_types": True })?
            raise ValueError(f"Multiple tables with different numbers of observations ({obs_type_to_num_rows[obs_type]} vs. {obs_num_rows}) are being used for obsType '{obs_type}'.")

        if feature_type not in feature_type_to_num_rows:
            feature_type_to_num_rows[feature_type] = var_num_rows
        elif feature_type_to_num_rows[feature_type] != var_num_rows:
            # TODO: same as above TODO.
            raise ValueError(f"Multiple tables with different numbers of variables ({feature_type_to_num_rows[feature_type]} vs. {var_num_rows}) are being used for featureType '{feature_type}'.")

        # TODO: configure all obsSets in the table here, to allow the user to select them regardless of the "color" parameter value,
        # rather than only when the "color" parameter is set to a categorical obs column (down below).
    else:
        # No annotating table exists, but if shapes, we can use the Parquet table to check the obs count.
        if element in sdata.shapes:
            obs_num_rows = sdata.shapes[element].shape[0]
            if obs_type not in obs_type_to_num_rows:
                obs_type_to_num_rows[obs_type] = obs_num_rows
            elif obs_type_to_num_rows[obs_type] != obs_num_rows:
                raise ValueError(f"Multiple tables with different numbers of observations ({obs_type_to_num_rows[obs_type]} vs. {obs_num_rows}) are being used for obsType '{obs_type}'.")

    obs_coordination = None
    feature_coordination = None
    if color is not None:
        if table_name is None:
            # color param must be a static color like "red" or "#FF0000"
            # TODO
            pass
        else:
            if color in sdata.tables[table_name].var.index: # gene
                feature_coordination = {
                    "obsType": obs_type,
                    "featureType": feature_type,
                    "featureSelection": [color],
                    "obsColorEncoding": "geneSelection",
                }

                if cmap is not None and cmap in ["viridis", "plasma", "jet", "greys"]:
                    feature_coordination["featureValueColormap"] = cmap
                elif cmap is None:
                    feature_coordination["featureValueColormap"] = "viridis"
                
                if norm is not None:
                    feature_coordination["featureValueColormapRange"] = [norm.vmin, norm.vmax]
            
            elif color in sdata.tables[table_name].obs.columns: # categorical?
                group_name = color.capitalize()

                # Configure the obsSets data wrapper properties.
                # Here we configure obsSets for wrapper_args
                wrapper_args["obs_set_paths"] = [f"tables/{table_name}/obs/{color}"]
                wrapper_args["obs_set_names"] = [group_name]

                obs_coordination = {
                    "obsType": obs_type,
                    "obsColorEncoding": "cellSetSelection",
                    "obsSetExpansion": [[group_name]],
                }

                # TODO: depends on https://github.com/vitessce/vitessce/issues/2254
                # obs_coordination["obsSetSelection"] = [[color]]
                if groups is not None:
                    obs_coordination["obsSetSelection"] = [
                        # Construct obs set paths.
                        [group_name, g] for g in groups
                    ]
                    if palette is not None:
                        if type(palette) is str:
                            # Broadcast single color to all groups.
                            palette = [palette for _ in groups]
                        elif type(palette) is list and len(groups) != len(palette):
                            raise ValueError("The length of 'groups' and 'palette' lists must be equal.")
                        
                        obs_coordination["obsSetColor"] = [
                            {
                                "path": [group_name, groups[i]],
                                "color": to_uint8_rgb(palette[i]),
                            } for i in range(len(groups))
                        ]
                else:
                    # Set to None, should initialize to children of first obsSet group by default.
                    obs_coordination["obsSetSelection"] = None
            else:
                # color param must be a static color like "red" or "#FF0000"
                # TODO
                pass

    return (obs_coordination, feature_coordination)


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
        self.shared_wrapper_args = {
            "sdata_path": self.sdata_filepath,
        }

        # TODO: Support same channel coordination across multiple layers, to render multiple image layers with linked channel settings?
        self.image_layers = [
            # Tuples of (wrapper_args, image_layer_coordination)
        ]
        self.segmentation_layers = [
            # Tuples of (wrapper_args, segmentation_layer_coordination, obs_coordination, feature_coordination)
        ]
        self.spot_layers = [
            # Tuples of (wrapper_args, spot_layer_coordination, obs_coordination, feature_coordination)
        ]
        self.point_layers = [
            # Tuples of (wrapper_args, point_layer_coordination)
        ]

        # For ensuring that counts of obs/var match if used for multiple layers.
        self.obs_type_to_num_rows = {}
        self.feature_type_to_num_rows = {}

    
    # References:
    # - https://spatialdata.scverse.org/projects/plot/en/latest/plotting.html#spatialdata_plot.pl.basic.PlotAccessor.render_images
    # - https://github.com/scverse/spatialdata-plot/blob/c9bae235c0521499fb4d1098b15c79619654e5dc/src/spatialdata_plot/pl/basic.py#L482
    def render_images(
            self,
            element=None,
            channel=None,
            cmap=None,
            norm=None,
            na_color=None,
            palette=None,
            alpha=1.0,
            **kwargs
        ):
        """
        Renders image.

        :param str element: location of image data inside "images" folder.
        :returns: Self, allows for chaining.
        """
        if not VitesscePlotAccessor._is_enabled:
            return self._pl.render_images(
                element=element,
                channel=channel,
                cmap=cmap,
                norm=norm,
                na_color=na_color,
                palette=palette,
                alpha=alpha,
                **kwargs,
            )

        # channel (list[str] | list[int] | str | int | None)
        #   To select specific channels to plot.
        #   Can be a single channel name/int or a list of channel names/ints.
        #   If None, all channels will be used.

        # cmap (list[Colormap | str] | Colormap | str | None)
        #   Colormap or list of colormaps for continuous annotations, see matplotlib.colors.Colormap.
        #   Each colormap applies to a corresponding channel.
 
        # palette (list[str] | str | None)
        #   Palette to color images.
        #   The number of palettes should be equal to the number of channels.

        # alpha (float | int, default 1.0)
        #   Alpha value for the images.
        #   Must be a numeric between 0 and 1.

        if type(channel) is list:
            # TODO: support lists of size 1 (broadcast/repeat to match num_channels length)?
            if type(cmap) is list:
                if len(channel) != len(cmap):
                    raise ValueError("The length of 'channel' and 'cmap' lists must be equal.")
            if type(palette) is list:
                if len(channel) != len(palette):
                    raise ValueError("The length of 'channel' and 'palette' lists must be equal.")
            if type(norm) is list:
                if len(channel) != len(norm):
                    raise ValueError("The length of 'channel' and 'norm' lists must be equal.")
            
        if element is None:
            # TODO: what does spatialdata-plot do in this case? use first image element? error if >1 images?
            raise ValueError("The 'element' parameter must be provided to render an image.")

        file_uid = f"image_{element}"
        wrapper_args = {
            "image_path": f"images/{element}",
            "coordination_values": {
                "fileUid": file_uid,
            }
        }

        # Palette logic in spatialdata-plot:
        # Reference: https://github.com/scverse/spatialdata-plot/blob/010560f7eebdd245693a8c55eede0f895a636f5c/src/spatialdata_plot/pl/utils.py#L685

        # RGB vs. non-RGB logic in spatialdata-plot:
        # Reference: https://github.com/scverse/spatialdata-plot/blob/010560f7eebdd245693a8c55eede0f895a636f5c/src/spatialdata_plot/pl/render.py#L865
        img = self.sdata.images[element]
        all_channels = img.coords["c"].values.tolist()
        img_dtype = img.dtype
        img_dtype_is_uint8 = img_dtype.kind == 'u' and img_dtype.itemsize == 1

        # Not ideal logic. Should ideally only use the OME-NGFF color model metadata. But this is what spatialdata-plot does.
        photometric_interpretation = "RGB" if palette is None and channel is None and len(all_channels) == 3 and img_dtype_is_uint8 else "BlackIsZero"

        # Configure image channel coordination.
        image_channel_coordination = None
        if channel is not None:
            # Normalize channels to a list.
            channel = [channel] if type(channel) in [int, str] else channel
            # Normalize palette to a list.
            if type(palette) is str:
                palette = [palette for _ in channel]
            if norm is not None and type(norm) is not list:
                norm = [norm for _ in channel]
            
            image_channel_coordination = [
                {
                    "spatialTargetC": ch,
                    **({ 'spatialChannelColor': to_uint8_rgb(palette[ch_i]) } if palette is not None else {}),
                    **({ 'spatialChannelWindow': [norm[ch_i].vmin, norm[ch_i].vmax] } if norm is not None else {}),
                    "spatialChannelVisible": True,
                }
                for ch_i, ch in enumerate(channel)
            ]
        
        # Configure image layer coordination.
        image_layer_coordination = {
            "fileUid": file_uid,
            "spatialLayerVisible": True,
            'spatialLayerOpacity': alpha,
            'photometricInterpretation': photometric_interpretation,
            **({
                'spatialLayerColormap': cmap
            } if cmap in ["viridis", "plasma", "jet", "greys"] else {}),
            **({} if na_color in [None, "default"] else {
                'spatialLayerTransparentColor': to_uint8_rgb(na_color)
            }),
            # Pass the image channel coordination if it was configured above.
            **({} if image_channel_coordination is None else {
                'imageChannel': image_channel_coordination
            }),
        }

        self.image_layers.append(
            (wrapper_args, image_layer_coordination)
        )

        return self.sdata
    
    # References:
    # - https://spatialdata.scverse.org/projects/plot/en/latest/plotting.html#spatialdata_plot.pl.basic.PlotAccessor.render_shapes
    # - https://github.com/scverse/spatialdata-plot/blob/c9bae235c0521499fb4d1098b15c79619654e5dc/src/spatialdata_plot/pl/basic.py#L156
    def render_shapes(self,
            element=None,
            color=None,
            fill_alpha=None,
            groups=None,
            palette=None,
            outline_width=None,
            outline_color=None,
            outline_alpha=None,
            cmap=None,
            norm=None,
            table_name=None,
            table_layer=None,
            **kwargs
        ):
        """
        Renders shapes, e.g. "cells".

        :param str element: location of shape data inside "shapes" folder.
        :param str color: gene.
        :param str cmap: color map (viridis, plasma, jet).
        :returns: Self, allows for chaining.
        """
        if not VitesscePlotAccessor._is_enabled:
            return self._pl.render_shapes(
                element=element,
                color=color,
                fill_alpha=fill_alpha,
                groups=groups,
                palette=palette,
                outline_width=outline_width,
                outline_color=outline_color,
                outline_alpha=outline_alpha,
                cmap=cmap,
                norm=norm,
                table_name=table_name,
                table_layer=table_layer,
                **kwargs
            )
        
        if element is None:
            # TODO: what does spatialdata-plot do in this case? use first shapes element? error if >1 shapes?
            raise ValueError("The 'element' parameter is required.")
        
        is_polygons = self.sdata.shapes[element]["geometry"].geom_type.iloc[0] == 'Polygon'

        file_uid = f"shapes_{element}"
        obs_type = "cell" if is_polygons else "spot"
        feature_type = "gene" # TODO: how to determine feature type? use heuristic based on num rows in table.var?

        wrapper_args = {
            "coordination_values": {
                "obsType": obs_type,
                "featureType": feature_type,
            }
        }

        # Vitessce only supports polygon and circle shapes.
        if is_polygons:
            wrapper_args["obs_segmentations_path"] = f"shapes/{element}"
            wrapper_args["coordination_values"]["fileUid"] = file_uid
            layer_coordination = {
                "fileUid": file_uid,
                "spatialLayerVisible": True,
                'segmentationChannel': [{
                    # We initialize with a single channel.
                    # SpatialData only supports single-channel segmentations.
                    "obsType": obs_type,
                    "featureType": feature_type,
                    "spatialChannelVisible": True,
                    "obsHighlight": None,
                }],
            }
        else:
            # Assume spots
            wrapper_args["obs_spots_path"] = f"shapes/{element}"
            layer_coordination = {
                "obsType": obs_type,
                "featureType": feature_type,
                "spatialLayerVisible": True,
                "obsHighlight": None,
            }

        # Shared coloring logic for polygons, spots, and labels.
        (obs_coordination, feature_coordination) = shared_render_shapes_and_labels(
            self.sdata, element, table_name, table_layer, color, cmap, norm, groups, palette, obs_type, feature_type,
            wrapper_args, self.obs_type_to_num_rows, self.feature_type_to_num_rows,
        )
        
        if is_polygons:
            self.segmentation_layers.append(
                (wrapper_args, layer_coordination, obs_coordination, feature_coordination)
            )
        else:
            self.spot_layers.append(
                (wrapper_args, layer_coordination, obs_coordination, feature_coordination)
            )
            if len(self.spot_layers) > 1:
                raise NotImplementedError("Multiple spot layers are not yet supported.")
        
        return self.sdata

    # References:
    # - https://spatialdata.scverse.org/projects/plot/en/latest/plotting.html#spatialdata_plot.pl.basic.PlotAccessor.render_labels
    # - https://github.com/scverse/spatialdata-plot/blob/c9bae235c0521499fb4d1098b15c79619654e5dc/src/spatialdata_plot/pl/basic.py#L598
    def render_labels(self,
            element=None,
            color=None,
            groups=None,
            palette=None,
            cmap=None,
            norm=None,
            outline_alpha=0.0,
            fill_alpha=0.4,
            table_name=None,
            table_layer=None,
            **kwargs
        ):
        """
        Renders label data.

        :param str element: location of label data in "labels" folder.
        :returns: Self, allows for chaining.
        """
        if not VitesscePlotAccessor._is_enabled:
            return self._pl.render_labels(
                element=element,
                color=color,
                groups=groups,
                palette=palette,
                cmap=cmap,
                norm=norm,
                outline_alpha=outline_alpha,
                fill_alpha=fill_alpha,
                table_name=table_name,
                table_layer=table_layer,
                **kwargs
            )

        if element is None:
            # TODO: what does spatialdata-plot do in this case? use first labels element? error if >1 labels?
            raise ValueError("The 'element' parameter must be provided to render labels.")
        
        file_uid = f"labels_{element}"
        obs_type = "cell"
        feature_type = "gene" # TODO: how to determine feature type? use heuristic based on num rows in table.var?

        wrapper_args = {
            "obs_segmentations_path": f"labels/{element}",
            "coordination_values": {
                "fileUid": file_uid,
                "obsType": obs_type,
                "featureType": feature_type,
            }
        }

        layer_coordination = {
            "fileUid": file_uid,
            "spatialLayerVisible": True,
            'segmentationChannel': [{
                # We initialize with a single channel.
                # SpatialData only supports single-channel segmentations.
                "obsType": obs_type,
                "featureType": feature_type,
                "spatialChannelVisible": True,
                "obsHighlight": None,
            }],
        }

        # Shared coloring logic for polygons, spots, and labels.
        (obs_coordination, feature_coordination) = shared_render_shapes_and_labels(
            self.sdata, element, table_name, table_layer, color, cmap, norm, groups, palette, obs_type, feature_type,
            wrapper_args, self.obs_type_to_num_rows, self.feature_type_to_num_rows,
        )
        
        self.segmentation_layers.append(
            (wrapper_args, layer_coordination, obs_coordination, feature_coordination)
        )
        
        return self.sdata

    # References:
    # - https://spatialdata.scverse.org/projects/plot/en/latest/plotting.html#spatialdata_plot.pl.basic.PlotAccessor.render_points
    # - https://github.com/scverse/spatialdata-plot/blob/c9bae235c0521499fb4d1098b15c79619654e5dc/src/spatialdata_plot/pl/basic.py#L338
    def render_points(self, element="", **kwargs):
        """
        Renders points.

        :param str element: location of point data in "points" folder.
        :returns: Self, allows for chaining.
        """
        if not VitesscePlotAccessor._is_enabled:
            return self._pl.render_points(element=element, **kwargs)

        file_uid = f"points_{element}"
        obs_type = "point"
        feature_type = "gene" # TODO: how to determine feature type? use heuristic based on num rows in table.var?

        wrapper_args = {
            "obs_points_path": f"points/{element}",
            "coordination_values": {
                "fileUid": file_uid,
                "obsType": obs_type,
                "featureType": feature_type,
            }
        }

        layer_coordination = {
            "obsType": obs_type,
            "obsHighlight": None,
            "fileUid": file_uid,
        }
        
        self.point_layers.append(
            (wrapper_args, layer_coordination)
        )

        return self.sdata
    
    def show(self, coordinate_systems=None, **kwargs):
        """
        Displays spatial plot.
        
        :returns: Vitessce widget. Learn more at the vitessce-python `docs <https://python-docs.vitessce.io/api_config.html#vitessce-widget>`_ .
        """
        if not VitesscePlotAccessor._is_enabled:
            return self._pl.show(coordinate_systems=coordinate_systems, **kwargs)
            
        self.vc = VitessceConfig(schema_version="1.0.18", name='SpatialData Plot')

        if not (coordinate_systems is None or isinstance(coordinate_systems, str)):
            raise NotImplementedError("A list of multiple 'coordinate_systems' is not yet supported.")

        dataset_uid = "A"
        dataset = self.vc.add_dataset(name='SpatialData Dataset', uid=dataset_uid)

        # TODO: de-duplicate wrapper_args if the same for multiple layers?
        for (layer_wrapper_args, _) in self.image_layers:
            img_wrapper = SpatialDataWrapper(**{
                **self.shared_wrapper_args,
                **({ "coordinate_system": coordinate_systems } if coordinate_systems is not None else {}),
                **layer_wrapper_args,
            })
            dataset = dataset.add_object(img_wrapper)
        
        for (layer_wrapper_args, _, _, _) in self.segmentation_layers:
            seg_wrapper = SpatialDataWrapper(**{
                **self.shared_wrapper_args,
                **({ "coordinate_system": coordinate_systems } if coordinate_systems is not None else {}),
                **layer_wrapper_args,
            })
            dataset = dataset.add_object(seg_wrapper)
        
        for (layer_wrapper_args, _, _, _) in self.spot_layers:
            spot_wrapper = SpatialDataWrapper(**{
                **self.shared_wrapper_args,
                **({ "coordinate_system": coordinate_systems } if coordinate_systems is not None else {}),
                **layer_wrapper_args,
            })
            dataset = dataset.add_object(spot_wrapper)
        
        for (layer_wrapper_args, _) in self.point_layers:
            points_wrapper = SpatialDataWrapper(**{
                **self.shared_wrapper_args,
                **({ "coordinate_system": coordinate_systems } if coordinate_systems is not None else {}),
                **layer_wrapper_args,
            })
            dataset = dataset.add_object(points_wrapper)


        # Add views (visualizations) to the configuration:
        spatial = self.vc.add_view("spatialBeta", dataset=dataset)
        layer_controller = self.vc.add_view("layerControllerBeta", dataset=dataset)
        obs_set_views = []
        feature_list_views = []

        obs_set_views_by_key = {}
        feature_list_views_by_key = {}


        # Collect all obs_coordination and feature_coordination information
        obs_coordination = []
        feature_coordination = []
        for (_, _, obs_coord, feature_coord) in self.segmentation_layers:
            if obs_coord is not None:
                obs_coordination.append(obs_coord)
            if feature_coord is not None:
                feature_coordination.append(feature_coord)
        for (_, _, obs_coord, feature_coord) in self.spot_layers:
            if obs_coord is not None:
                obs_coordination.append(obs_coord)
            if feature_coord is not None:
                feature_coordination.append(feature_coord)

        # Add obsSet and featureList views.
        for obs_coord in obs_coordination:
            obs_set_view = self.vc.add_view("obsSets", dataset=dataset)
            obs_set_views.append(obs_set_view)
            obs_set_views_by_key[obs_coord["obsType"]] = obs_set_view
        for feature_coord in feature_coordination:
            feature_list_view = self.vc.add_view("featureList", dataset=dataset)
            feature_list_views.append(feature_list_view)
            feature_list_views_by_key[feature_coord["featureType"]] = feature_list_view

        spatial_views = [spatial, layer_controller]
        control_views = [layer_controller, *obs_set_views, *feature_list_views]
        all_views = [spatial, *control_views]

        # Coordinate views.
        self.vc.link_views_by_dict(spatial_views, {
            "imageLayer": CL([
                {
                    **layer_dict,
                    **({} if "imageChannel" not in layer_dict else {
                        'imageChannel': CL([
                            {
                                **channel_dict,
                            }
                            for channel_dict in layer_dict['imageChannel']
                        ])
                    })
                }
                for (_, layer_dict) in self.image_layers
            ]),
        }, meta=True, scope_prefix=get_initial_coordination_scope_prefix(dataset_uid, "image"))

        # Collect per-obsType coordination types and scopes.
        obs_coordination_by_key = {}
        for obs_coord in obs_coordination:
            coordination_key = obs_coord["obsType"] # TODO: is this the best key to use?

            # Create coordination scope objects.
            ct_names = []
            ct_vals = []
            for ct_name, ct_val in obs_coord.items():
                ct_names.append(ct_name)
                ct_vals.append(ct_val)
            
            ct_scopes = self.vc.add_coordination(*ct_names)
            for i, ct_scope in enumerate(ct_scopes):
                ct_scope.set_value(ct_vals[i])
            obs_coordination_by_key[coordination_key] = dict(zip(ct_names, ct_scopes))
        
        feature_coordination_by_key = {}
        for feature_coord in feature_coordination:
            coordination_key = feature_coord["featureType"] # TODO: is this the best key to use?

            # Create coordination scope objects.
            ct_names = []
            ct_vals = []
            for ct_name, ct_val in feature_coord.items():
                ct_names.append(ct_name)
                ct_vals.append(ct_val)
            
            ct_scopes = self.vc.add_coordination(*ct_names)
            for i, ct_scope in enumerate(ct_scopes):
                ct_scope.set_value(ct_vals[i])
            feature_coordination_by_key[coordination_key] = dict(zip(ct_names, ct_scopes))
        
        self.vc.link_views_by_dict(spatial_views, {
            "segmentationLayer": CL([
                {
                    **layer_dict,
                    'segmentationChannel': CL([
                        {
                            **channel_dict,
                            **obs_coordination_by_key.get(channel_dict.get("obsType"), {}),
                            **feature_coordination_by_key.get(channel_dict.get("featureType"), {}),
                        }
                        for channel_dict in layer_dict.get('segmentationChannel', [{}])
                    ])
                }
                for (_, layer_dict, _, _) in self.segmentation_layers
            ]),
        }, meta=True, scope_prefix=get_initial_coordination_scope_prefix(dataset_uid, "obsSegmentations"))

        self.vc.link_views_by_dict(spatial_views, {
            "spotLayer": CL([
                # TODO: should this be restricted to a single dict (not array), since Vitessce supports only a single spot layer?
                {
                    **layer_dict,
                    **obs_coordination_by_key.get(layer_dict.get("obsType"), {}),
                    **feature_coordination_by_key.get(layer_dict.get("featureType"), {}),
                }
                for (_, layer_dict, _, _) in self.spot_layers
            ]),
        }, meta=True, scope_prefix=get_initial_coordination_scope_prefix(dataset_uid, "obsSpots"))

        self.vc.link_views_by_dict(spatial_views, {
            "pointLayer": CL([
                {
                    **layer_dict,
                    **obs_coordination_by_key.get(layer_dict.get("obsType"), {}),
                    **feature_coordination_by_key.get(layer_dict.get("featureType"), {}),
                }
                for (_, layer_dict) in self.point_layers
            ]),
        }, meta=True, scope_prefix=get_initial_coordination_scope_prefix(dataset_uid, "obsPoints"))

        # Set up coordination for control views.
        for key, obs_set_view in obs_set_views_by_key.items():
            self.vc.link_views_by_dict([obs_set_view], obs_coordination_by_key.get(key, {}), meta=False)

        for key, feature_list_view in feature_list_views_by_key.items():
            self.vc.link_views_by_dict([feature_list_view], feature_coordination_by_key.get(key, {}), meta=False)
        
        
        # Layout the views
        self.vc.layout(hconcat(spatial, vconcat(*control_views), split=[2, 1]))
        
        vw = _to_widget(self.vc)

        # Cleanup
        self._init_params()

        return vw
