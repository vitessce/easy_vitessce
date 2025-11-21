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
from easy_vitessce.data import _get_sdata_wrapper_params

# Internal function for shared logic between render_shapes and render_labels.
def _shared_render_shapes_and_labels(
        sdata, element, table_name, table_layer, color, cmap, norm, groups, palette, obs_type, feature_type, is_spots, fill_alpha, outline_alpha, outline_width, outline_color,
        # Note: These dict params are modified by this function.
        wrapper_args, obs_type_to_num_rows, feature_type_to_num_rows,
    ):

    extra_layer_coordination = {}

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
        # Update: commented out since it causes issues as obs_type is always "cell" but the element is different.
        """
        if element in sdata.shapes:
            obs_num_rows = sdata.shapes[element].shape[0]
            if obs_type not in obs_type_to_num_rows:
                obs_type_to_num_rows[obs_type] = obs_num_rows
            elif obs_type_to_num_rows[obs_type] != obs_num_rows:
                raise ValueError(f"Multiple tables with different numbers of observations ({obs_type_to_num_rows[obs_type]} vs. {obs_num_rows}) are being used for obsType '{obs_type}'.")
        """
        pass

    obs_coordination = None
    feature_coordination = None
    is_maybe_static_color = False
    if color is not None:
        if table_name is None:
            # color param must be a static color like "red" or "#FF0000"
            is_maybe_static_color = True
        else:
            if color in sdata.tables[table_name].var.index: # gene
                extra_layer_coordination["obsColorEncoding"] = "geneSelection"
                feature_coordination = {
                    "obsType": obs_type,
                    "featureType": feature_type,
                    "featureSelection": [color],
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

                extra_layer_coordination["obsColorEncoding"] = "cellSetSelection"
                obs_coordination = {
                    "obsType": obs_type,
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
                is_maybe_static_color = True
    else:
        is_maybe_static_color = True
    
    if is_maybe_static_color:
        extra_layer_coordination["obsColorEncoding"] = "spatialChannelColor" if not is_spots else "spatialLayerColor"
        if color is not None:
            extra_layer_coordination["spatialChannelColor" if not is_spots else "spatialLayerColor"] = to_uint8_rgb(color)
    
    # Handling of alpha/fill settings.
    # We can only support these things partially.
    # - Case 1: fill_alpha == outline_alpha
    # - Case 2: fill_alpha != outline_alpha but fill_alpha is 0.0 (effectively stroked)
    # - Case 3: fill_alpha != outline_alpha and both > 0.0 (requires both fill and stroke colors with different alphas)
    # Case 3 is not currently supported by Vitessce unless we create two separate layers (one for fill, one for outline).
    if fill_alpha == 0.0:
        # Stroked.
        if is_spots:
            extra_layer_coordination["spatialSpotFilled"] = False
        else:
            extra_layer_coordination["spatialSegmentationFilled"] = False
    
    if outline_width is not None:
        if is_spots:
            extra_layer_coordination["spatialSpotStrokeWidth"] = outline_width
        else:
            extra_layer_coordination["spatialSegmentationStrokeWidth"] = outline_width

    layer_alpha = None
    if outline_alpha is not None and (fill_alpha == outline_alpha or fill_alpha == 0.0):
        layer_alpha = outline_alpha
    elif fill_alpha is not None:
        layer_alpha = fill_alpha

    if layer_alpha is not None:
        if is_spots:
            extra_layer_coordination["spatialLayerOpacity"] = layer_alpha
        else:
            extra_layer_coordination["spatialChannelOpacity"] = layer_alpha

    return (extra_layer_coordination, obs_coordination, feature_coordination)


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

        # This is the static PlotAccessor instance that will be used when monkeypatching is not enabled.
        self._pl = PlotAccessor(sdata)

        self.did_init = False
        self._maybe_init()
    
    def _maybe_init(self):
        # We cannot assume that __init__ has been called as expected,
        # for instance if ev.enable_plots is called after creating the SpatialData object.
        # Instead, we call _maybe_init at the start of each public method.

        if not self.did_init:
            self._init_params()
            self.did_init = True

    def _init_params(self):
        # Initialize or re-initialize plotting state.
        # Called in constructor (to initialize) and after .pl.show() (to clean up state prior to next .pl.render_something).
        self.shared_wrapper_args = _get_sdata_wrapper_params(self.sdata)

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

        :param str element: Name of the image element.
        :param channel: To select specific channels to plot.
        :type channel: list[str] or list[int] or str or int or None
        :param cmap: Colormap name such as "viridis".
        :type cmap: str or None
        :param norm: Normalization or list of normalizations for continuous annotations.
        :type norm: list[matplotlib.colors.Normalize] or matplotlib.colors.Normalize or None
        :param na_color: Color that should be rendered as transparent.
        :type na_color: str or None
        :param palette: Palette to color images. If list, the number of colors should be equal to the number of channels.
        :type palette: list[str] or str or None
        :param alpha: Alpha value for the images, between 0.0 and 1.0. By default, 1.0.
        :type alpha: float or int
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
        
        self._maybe_init()

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
        try:
            all_channels = img.coords["c"].values.tolist()
        except KeyError:
            # TODO: use a better way than try/except of determining whether this is a multi-resolution image.
            all_channels = img.scale0.coords["c"].values.tolist()
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

        :param str element: Name of the shapes element.
        :param color: Name of an obs column, var index value, or color-like string.
        :type color: str or None
        :param fill_alpha: Alpha value for filling shapes, between 0.0 and 1.0.
        :type fill_alpha: float or int or None
        :param groups: List of obs group names to select.
        :type groups: list[str] or str or None
        :param palette: Palette to color shapes. If list, the number of colors should be equal to the number of groups.
        :type palette: list[str] or str or None
        :param outline_width: Width of the shape outlines.
        :type outline_width: float or int or None
        :param outline_alpha: Alpha value for shape outlines, between 0.0 and 1.0.
        :type outline_alpha: float or int or None
        :param cmap: Quantitative colormap name, such as "viridis".
        :type cmap: str or None
        :param norm: Normalization for quantitative colormap.
        :type norm: matplotlib.colors.Normalize or None
        :param table_name: Name of an annotating table to use for coloring.
        :type table_name: str or None
        :param table_layer: Name of the layer in the annotating table to use for coloring.
        :type table_layer: str or None
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
        
        self._maybe_init()
        
        if element is None:
            # TODO: what does spatialdata-plot do in this case? use first shapes element? error if >1 shapes?
            raise ValueError("The 'element' parameter is required.")
        
        is_polygons = self.sdata.shapes[element]["geometry"].geom_type.iloc[0] == 'Polygon'
        is_spots = not is_polygons

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
        (extra_layer_coordination, obs_coordination, feature_coordination) = _shared_render_shapes_and_labels(
            self.sdata, element, table_name, table_layer, color, cmap, norm, groups, palette, obs_type, feature_type, is_spots, fill_alpha, outline_alpha, outline_width, outline_color,
            wrapper_args, self.obs_type_to_num_rows, self.feature_type_to_num_rows,
        )
        
        if is_polygons:
            layer_coordination["segmentationChannel"][0].update(extra_layer_coordination)
            self.segmentation_layers.append(
                (wrapper_args, layer_coordination, obs_coordination, feature_coordination)
            )
        else:
            layer_coordination.update(extra_layer_coordination)
            self.spot_layers.append(
                (wrapper_args, layer_coordination, obs_coordination, feature_coordination)
            )
        
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

        :param str element: Name of the labels element.
        :param color: Name of an obs column, var index value, or color-like string.
        :type color: str or None
        :param groups: List of obs group names to select.
        :type groups: list[str] or str or None
        :param palette: Palette to color labels. If list, the number of colors should be equal to the number of groups.
        :type palette: list[str] or str or None
        :param cmap: Quantitative colormap name, such as "viridis".
        :type cmap: str or None
        :param norm: Normalization for quantitative colormap.
        :type norm: matplotlib.colors.Normalize or None
        :param outline_alpha: Alpha value for label outlines, between 0.0 and 1.0.
        :type outline_alpha: float or int or None
        :param fill_alpha: Alpha value for filling labels, between 0.0 and 1.0.
        :type fill_alpha: float or int or None
        :param table_name: Name of an annotating table to use for coloring.
        :type table_name: str or None
        :param table_layer: Name of the layer in the annotating table to use for coloring.
        :type table_layer: str or None
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
        
        self._maybe_init()

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

        is_spots = False
        outline_width = None # Not supported for labels in spatialdata-plot.
        outline_color = None # Not supported for labels in spatialdata-plot.

        # Shared coloring logic for polygons, spots, and labels.
        (extra_layer_coordination, obs_coordination, feature_coordination) = _shared_render_shapes_and_labels(
            self.sdata, element, table_name, table_layer, color, cmap, norm, groups, palette, obs_type, feature_type, is_spots, fill_alpha, outline_alpha, outline_width, outline_color,
            wrapper_args, self.obs_type_to_num_rows, self.feature_type_to_num_rows,
        )

        layer_coordination["segmentationChannel"][0].update(extra_layer_coordination)
        
        self.segmentation_layers.append(
            (wrapper_args, layer_coordination, obs_coordination, feature_coordination)
        )
        
        return self.sdata

    # References:
    # - https://spatialdata.scverse.org/projects/plot/en/latest/plotting.html#spatialdata_plot.pl.basic.PlotAccessor.render_points
    # - https://github.com/scverse/spatialdata-plot/blob/c9bae235c0521499fb4d1098b15c79619654e5dc/src/spatialdata_plot/pl/basic.py#L338
    def render_points(self, element=None, **kwargs):
        """
        Renders points.

        :param str element: Name of points element.
        :returns: Self, allows for chaining.
        """
        if not VitesscePlotAccessor._is_enabled:
            return self._pl.render_points(element=element, **kwargs)
        
        self._maybe_init()

        if element is None:
            # TODO: what does spatialdata-plot do in this case? use first points element? error if >1 points?
            raise ValueError("The 'element' parameter must be provided to render points.")

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
        if len(self.image_layers) > 0:
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
        
        if len(self.segmentation_layers) > 0:
            self.vc.link_views_by_dict(spatial_views, {
                "segmentationLayer": CL([
                    {
                        **layer_dict,
                        'segmentationChannel': CL([
                            {
                                **feature_coordination_by_key.get(channel_dict.get("featureType"), {}),
                                **obs_coordination_by_key.get(channel_dict.get("obsType"), {}),
                                **channel_dict,
                            }
                            for channel_dict in layer_dict.get('segmentationChannel', [{}])
                        ])
                    }
                    for (_, layer_dict, _, _) in self.segmentation_layers
                ]),
            }, meta=True, scope_prefix=get_initial_coordination_scope_prefix(dataset_uid, "obsSegmentations"))

        if len(self.spot_layers) > 0:
            self.vc.link_views_by_dict(spatial_views, {
                "spotLayer": CL([
                    {
                        **feature_coordination_by_key.get(layer_dict.get("featureType"), {}),
                        **obs_coordination_by_key.get(layer_dict.get("obsType"), {}),
                        **layer_dict,
                    }
                    for (_, layer_dict, _, _) in self.spot_layers
                ]),
            }, meta=True, scope_prefix=get_initial_coordination_scope_prefix(dataset_uid, "obsSpots"))

        if len(self.point_layers) > 0:
            self.vc.link_views_by_dict(spatial_views, {
                "pointLayer": CL([
                    {
                        **feature_coordination_by_key.get(layer_dict.get("featureType"), {}),
                        **obs_coordination_by_key.get(layer_dict.get("obsType"), {}),
                        **layer_dict,
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
