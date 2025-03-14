from importlib.resources import files


def get_example_path(file_name):
    """Get full path to given example file or directory.

    :param file_name: Name of the example file or directory
    :type file_name: str
    :return: Path to example file or directory
    :rtype: str

    :Example:

    >>> import tcrconvert
    >>> tcrconvert.get_example_path('tenx.csv') # doctest: +ELLIPSIS
    '.../tcrconvert/examples/tenx.csv'
    """

    out = files('tcrconvert').joinpath('examples', file_name)
    return str(out)
