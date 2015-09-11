#!/bin/bash

# ----------------------------------------------------------------------------
# If not defined... define the base environment variables

# Where the software will be installed
if [ -z "$PREFIX" ]; then
    export PREFIX=/usr/local
fi

# Where the external libraries are installed
if [ -z "$COTS" ]; then
    export COTS=/data/cots
fi

# ----------------------------------------------------------------------------
# Always required to be installed for building since it is part of the ESPA
export ESPAINC=${PREFIX}/espa-common/include
export ESPALIB=${PREFIX}/espa-common/lib

# Setup base paths to external libraries
libxml2_path=${COTS}

export XML2INC=${libxml2_path}/libxml2/include/libxml2
export XML2LIB=${libxml2_path}/libxml2/lib

