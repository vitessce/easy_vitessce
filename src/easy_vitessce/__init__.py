from .configure_plots import configure_plots, enable_plots, disable_plots
from .data import register_data_params, register_data_path, register_data_store, register_data_url
from .widget import config

# Side effect: Enable plotting functions on import.
enable_plots()