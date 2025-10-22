from matplotlib.colors import to_rgb

def to_uint8_rgb(color_str):
    """
    Converts a color string to a list of RGB uint8 values.

    :param str color_str: The color string, e.g., "red" or "#ff0000".
    :return: A list of three integers [R, G, B] from 0 to 255.
    :rtype: list[int]
    """
    rgb_float = to_rgb(color_str)
    return [int(c * 255) for c in rgb_float]