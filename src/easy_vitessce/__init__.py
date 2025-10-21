from .configure_plots import configure_plots, enable_plots, disable_plots
from .data import register_data_path
from .widget import config

# Side effect: Enable plotting functions on import.
enable_plots()