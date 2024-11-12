__version__ = '26.0.0'


def is_augur_version_compatible(version):
    """
    Checks if the provided **version** is the same major version
    as the currently running version of augur.

    Parameters
    ----------
    version : str
        version to check against the current version

    Returns
    -------
    Bool

    """
    import packaging.version

    current_version = packaging.version.parse(__version__)
    this_version = packaging.version.parse(version)
    return this_version.release[0] == current_version.release[0]
