
# Configuration for logging in to the UCAR system
ucar_login_credentials = {
    'login_url': 'https://rda.ucar.edu/cgi-bin/login',
    'login_data': {
        'email': None,   # 'xxx@yyy.zzz'
        'passwd': None,  # 'XXXXXXXXX'
        'action': 'login',
    }
}

# URL location for the NARR data on the UCAR system
# 0: year
# 1: filename
UCAR_URL_FORMAT = 'http://rda.ucar.edu/data/ds608.0/3HRLY/{0}/{1}'

# The format describing the name of the files on the remote UCAR system
# 0: year
# 1: month
# 2: start day (inclusive)
# 3: end day (inclusive)
REMOTE_NAME_FORMAT = 'NARR3D_{0:04}{1:02}_{2:02}{3:02}'

# Our internal naming scheme
# 0: variable (HGT, SPFH, TMP)
# 1: year
# 2: month
# 3: day
# 4: hour (0000, 0300, ..., 2100)
# 5: file extension (hdr, grb)
ARCHIVE_NAME_FORMAT = 'NARR_3D.{0}.{1:04}{2:02}{3:02}.{4:04}.{5}'

# Our internal location
# 0: base auxillary directory
# 1: year
# 2: month
# 3: day
ARCHIVE_DIRECTORY_FORMAT = '{0}/{1:0>4}/{2:0>2}/{3:0>2}'

# The NARR variable this application requires
NARR_VARIABLES = ['HGT', 'TMP', 'SPFH']
