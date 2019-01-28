"""Provide configuration file reading
"""


import os
from configparser import ConfigParser
from collections import namedtuple

from st_aux_exception import AuxiliaryError


class AuxiliaryConfigError(AuxiliaryError):
    """Config errors
    """
    pass


CONFIG_FILENAME = 'st_merra_auxiliary.conf'
CONFIG_DIR = '.usgs/espa'


def config_file_full_path():
    """Build the full path to the config file

    Raises:
        AuxiliaryConfigError(message)
    """

    # Use the user's home directory as the base source directory for
    # configuration
    if 'HOME' not in os.environ:
        raise AuxiliaryConfigError('[HOME] not found in environment')
    home_dir = os.environ.get('HOME')

    # Build the full path to the configuration file
    config_path = os.path.join(home_dir, CONFIG_DIR, CONFIG_FILENAME)

    return config_path


def read_config():
    """Reads a config file

    Raises:
        AuxiliaryConfigError(message)
    """

    config_path = config_file_full_path()

    if not os.path.isfile(config_path):
        raise AuxiliaryConfigError('Missing configuration file [{}]'
                                    .format(config_path))

    # Create the object and load the configuration
    cfg = ConfigParser()
    cfg.read(config_path)

    return cfg


NASAInfo = namedtuple('NASAInfo', ('data_url_format',
                                   'data_name_format',
                                   'login_url',
                                   'action',
                                   'username',
                                   'password'))


ConfigInfo = namedtuple('ConfigInfo', ('archive_directory_format',
                                       'archive_name_format',
                                       'base_archive_directory',
                                       'search_date_range',
                                       'transfer_block_size',
                                       'nasa'))


def get_config():
    """Gets the configuration and converts it

    Raises:
        AuxiliaryConfigError(message)
    """

    cfg = read_config()

    merra_section = 'merra2'
    nasa_section = 'nasa'

    nasa = NASAInfo(data_url_format=cfg.get(nasa_section, 'data_url_format'),
                    data_name_format=cfg.get(nasa_section, 'data_name_format'),
                    login_url=cfg.get(nasa_section, 'login_url'),
                    action=cfg.get(nasa_section, 'action'),
                    username=cfg.get(nasa_section, 'username'),
                    password=cfg.get(nasa_section, 'password'))

    return ConfigInfo(archive_directory_format=cfg.get(merra_section, 'archive_directory_format'),
                      archive_name_format=cfg.get(merra_section, 'archive_name_format'),
                      base_archive_directory=cfg.get(merra_section, 'base_archive_directory'),
                      search_date_range=cfg.get(merra_section, 'search_date_range'),
                      transfer_block_size=int(cfg.get(merra_section, 'transfer_block_size')),
                      nasa=nasa)
