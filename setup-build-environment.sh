#!/bin/bash

# Always required to be installed for building
export ESPAINC=${PREFIX}/espa-common/include
export ESPALIB=${PREFIX}/espa-common/lib

# Setup base paths to external libraries
libxml2_path=${COTS}

export XML2INC=${libxml2_path}/libxml2/include/libxml2
export XML2LIB=${libxml2_path}/libxml2/lib

