from donfig import Config

config = Config('easy_vitessce', defaults=[{
    'widget': {
        # Parameters of vc.widget().
        'js_dev_mode': True,
    }
}])

def to_widget(vc):
    """
    Converts a VitessceConfig object to a widget.

    :param vc: A VitessceConfig object.
    :returns: A widget representation of the configuration.
    """
    return vc.widget(
        **config.get('widget'),
    )