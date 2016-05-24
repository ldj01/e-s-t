
'''
    File: lst_build_points.py

    Purpose: Builds a directory structure of points and required information
             to be used for generation and interpolation of atmospheric
             information by follow on applications.

    Project: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    License: NASA Open Source Agreement 1.3
'''


class LSTError(Exception):
    """Base exception for LST"""
    pass


class MissingBandError(LSTError):
    """Exception to use for missing bands"""
    pass
