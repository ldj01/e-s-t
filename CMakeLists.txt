cmake_minimum_required(VERSION 3.6)

project(espa-surface-temperature)

# Set the package version.
include(package_version.cmake)

include(FindPkgConfig)

if (NOT XML_INCLUDES)
    find_package(LibXml2)
    if (LIBXML2_FOUND)
        set (XML_INCLUDES ${LIBXML2_INCLUDE_DIR})
        set (XML_LIBS ${LIBXML2_LIBRARIES})
    endif ()
endif ()
set (XML_INCLUDES ${XML_INCLUDES} CACHE PATH "XML2 header file path" FORCE)
set (XML_LIBS ${XML_LIBS} CACHE STRING "XML2 libraries" FORCE)

# Check whether the HDFEOS flags or environment variables are set. If so,
# then use the HDFEOS GCTP library.  Otherwise, use the GCTP2 library.
if (NOT HDFEOS_LIBS)
    if (DEFINED ENV{HDFEOS_LIB})
        set (GCTP2_LIBS -L$ENV{HDFEOS_LIB} -lGctp -lm)
        set (USE_HDFEOS TRUE)
    elseif (USE_HDFEOS)
        pkg_check_modules(HDFEOS hdf-eos4)
        if (HDFEOS_FOUND)
            set (HDFEOS_LIBS ${HDFEOS_LDFLAGS})
        endif ()
    endif ()
else ()
    set (USE_HDFEOS TRUE)
endif ()
if (USE_HDFEOS)
    message(STATUS "Using HDFEOS library.")
    if (DEFINED ENV{HDFEOS_LIB})
        set (GCTP2_LIBS "-L$ENV{HDFEOS_LIB} -lGctp -lm")
    else ()
        set (GCTP2_LIBS "${HDFEOS_LIBS} -lGctp -lm")
    endif ()
else ()
    message(STATUS "Using GCTP2 package.")
    if (NOT GCTP2_LIBS)
        set (GCTP2_LIBS "-L$ENV{GCTP2LIB} -lgctp -lm")
    endif ()
endif ()
set (GCTP2_LIBS ${GCTP2_LIBS} CACHE STRING "GCTP2 libraries" FORCE)

# Set the paths for the ESPA headers and libraries.
if (NOT ESPA_INCLUDES)
    if (DEFINED ENV{ESPAINC})
        set (ESPA_INCLUDES $ENV{ESPAINC})
    else ()
        set (ESPA_INCLUDES ${CMAKE_INSTALL_PREFIX}/include)
    endif ()
endif ()
if (NOT ESPA_LIBS)
    if (DEFINED ENV{ESPALIB})
        set (ESPA_LIBS $ENV{ESPALIB})
    else ()
        set (ESPA_LIBS ${CMAKE_INSTALL_PREFIX}/lib)
    endif ()
endif ()
set (ESPA_INCLUDES ${ESPA_INCLUDES} CACHE PATH "ESPA header file path" FORCE)
set (ESPA_LIBS ${ESPA_LIBS} CACHE PATH "ESPA library path" FORCE)

enable_testing()

# Create a custom target for the RIT auxiliary data tool installation.
add_custom_target(rit-aux
                  COMMAND "${CMAKE_COMMAND}"
                          -DCMAKE_INSTALL_COMPONENT=rit-aux
                          -P "${CMAKE_BINARY_DIR}/cmake_install.cmake")

# Create "make install-data" target.
add_custom_target(install-data
                  COMMAND "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=data
                          -P "${CMAKE_BINARY_DIR}/cmake_install.cmake")

add_subdirectory(scripts)
add_subdirectory(provisional-st)
add_subdirectory(st_auxiliary_data)
