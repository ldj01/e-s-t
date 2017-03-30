"""Provide configuration file reading
"""


import os
from ConfigParser import ConfigParser
from collections import namedtuple

from lst_aux_exception import AuxiliaryError


class AuxiliaryConfigError(AuxiliaryError):
    """Config errors
    """
    pass


CONFIG_FILENAME = 'lst_auxiliary.conf'
CONFIG_DIR = '.usgs/espa'


def config_file_full_path():
    """Build the full path to the config file

    Raises:
        AuxiliaryConfigError(message)
    """

    # Use the users home directory as the base source directory for
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


UCARInfo = namedtuple('UCARInfo', ('data_url_format',
                                   'data_name_format',
                                   'login_url',
                                   'action',
                                   'email',
                                   'passwd'))


NCEPInfo = namedtuple('NCEPInfo', ('data_url_format',
                                   'data_name_format'))


ConfigInfo = namedtuple('ConfigInfo', ('archive_directory_format',
                                       'archive_name_format',
                                       'base_archive_directory',
                                       'search_date_range',
                                       'transfer_block_size',
                                       'ucar',
                                       'ncep'))


def get_config():
    """Gets the configuration and converts it

    Raises:
        AuxiliaryConfigError(message)
    """

    cfg = read_config()

    narr_section = 'narr'
    ucar_section = 'ucar'
    ncep_section = 'ncep'

    ucar = UCARInfo(data_url_format=cfg.get(ucar_section, 'data_url_format'),
                    data_name_format=cfg.get(ucar_section, 'data_name_format'),
                    login_url=cfg.get(ucar_section, 'login_url'),
                    action=cfg.get(ucar_section, 'action'),
                    email=cfg.get(ucar_section, 'email'),
                    passwd=cfg.get(ucar_section, 'passwd'))

    ncep = NCEPInfo(data_url_format=cfg.get(ncep_section, 'data_url_format'),
                    data_name_format=cfg.get(ncep_section, 'data_name_format'))

    return ConfigInfo(archive_directory_format=cfg.get(narr_section, 'archive_directory_format'),
                      archive_name_format=cfg.get(narr_section, 'archive_name_format'),
                      base_archive_directory=cfg.get(narr_section, 'base_archive_directory'),
                      search_date_range=cfg.get(narr_section, 'search_date_range'),
                      transfer_block_size=int(cfg.get(narr_section, 'transfer_block_size')),
                      ucar=ucar,
                      ncep=ncep)
