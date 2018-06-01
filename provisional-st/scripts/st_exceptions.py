
'''
    File: st_build_points.py

    Purpose: Builds a directory structure of points and required information
             to be used for generation and interpolation of atmospheric
             information by follow on applications.

    Project: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    License: NASA Open Source Agreement 1.3
'''


class STError(Exception):
    """Base exception for ST"""
    pass


class MissingBandError(STError):
    """Exception to use for missing bands"""
    pass

class NoTilesError(STError):
    """Exception to use for not finding any tiles"""
    pass

class InaccessibleTileError(STError):
    """Exception to use for not being able to access a tile"""
    pass
