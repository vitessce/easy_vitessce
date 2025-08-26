from donfig import Config

config = Config('easy_vitessce', defaults=[{
    'widget_function': 'widget', # Options: 'widget' or 'display'
    'widget': {
        # Parameters of vc.widget().
        'js_dev_mode': True,
    },
    'data': {
        # Configure the directory used for saving local anndata/spatialdata objects on-disk.
        'out_dir': 'data',
        'overwrite': False,
        'anndata_format': 'zarr',  # Options: 'zarr' or 'h5ad'
        'wrapper_param_suffix': '_path' # Options: '_path' or '_store' or '_url'
    },

    # Per-plot configurations that are Vitessce-specific?
    # For example, for `embedding`, whether to include a featureList view when `color` is a cell set. Or whether to include an obsSets view when `color` is a gene.

}])

def _to_widget(vc):
    """
    Converts a VitessceConfig object to a widget.

    :param vc: A VitessceConfig object.
    :returns: A widget representation of the configuration.
    """
    widget_params = config.get('widget')
    if config.get('widget_function') == 'display':
        return vc.display(**widget_params)
    else:
        return vc.widget(**widget_params)