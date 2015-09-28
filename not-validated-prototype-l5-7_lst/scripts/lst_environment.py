
'''
    FILE: lst_environment.py

    PURPOSE: Validates and provides access to environment variables.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3

    HISTORY:

    Date              Reason
    ----------------  --------------------------------------------------------
    Sept/2015          Initial implementation
'''


import os


# ----------------------------------------------------------------------------
class Environment(object):
    '''
    Description:
        Grab environment variables from the system which are used during LST
        processing.  Some of them are very simply validated.
    '''

    def __init__(self):
        super(Environment, self).__init__()

        self.lst_data_dir = os.environ.get('LST_DATA_DIR')
        if self.lst_data_dir is None:
            raise Exception('Missing environment variable LST_DATA_DIR')

        # Verify that the lst_data_dir directory exists
        if not os.path.isdir(self.lst_data_dir):
            raise Exception('LST_DATA_DIR directory does not exist')

        self.lst_aux_dir = os.environ.get('LST_AUX_DIR')
        if self.lst_aux_dir is None:
            raise Exception('Missing environment variable LST_AUX_DIR')

        # Verify that the lst_aux_dir directory exists
        if not os.path.isdir(self.lst_aux_dir):
            raise Exception('LST_AUX_DIR directory does not exist')

        self.aster_ged_server_name = os.environ.get('ASTER_GED_SERVER_NAME')
        if self.aster_ged_server_name is None:
            raise Exception('Missing environment variable'
                            ' ASTER_GED_SERVER_NAME')

    def get_lst_data_directory(self):
        '''Returns the LST data directory'''
        return self.lst_data_dir

    def get_lst_aux_directory(self):
        '''Returns the LST data directory'''
        return self.lst_aux_dir

    def get_aster_ged_server_name(self):
        '''Returns the ASTER GED server name'''
        return self.aster_ged_server_name
