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
# Always required to be installed for building since it is part of almost
# every ESPA science application
if [ -z "$ESPAINC" ]; then
    export ESPAINC=${PREFIX}/espa-common/include
    export ESPALIB=${PREFIX}/espa-common/lib
fi

# Setup base paths to external libraries
base_path=${COTS}

if [ -z "$XML2INC" ]; then
    export XML2INC=${base_path}/libxml2/include/libxml2
    export XML2LIB=${base_path}/libxml2/lib
fi

if [ -z "$LZMALIB" ]; then
    export LZMAINC=${base_path}/xz/include
    export LZMALIB=${base_path}/xz/lib
fi

if [ -z "$ZLIBLIB" ]; then
    export ZLIBINC=${base_path}/zlib/include
    export ZLIBLIB=${base_path}/zlib/lib
fi

