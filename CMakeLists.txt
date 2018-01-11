cmake_minimum_required(VERSION 3.6)

project(espa-surface-temperature)
set (${PROJECT_NAME}_VERSION_MAJOR 1)
set (${PROJECT_NAME}_VERSION_MINOR 0)
set (${PROJECT_NAME}_VERSION_PATCH 0)
set (VERSION ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH})

set (PACKAGE ${PROJECT_NAME})

set(CPACK_PACKAGE_VERSION ${VERSION})
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "/build/;/.git/;/.gz/;~$;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)

if (XML2INC)
    set (XML_INCLUDES ${XML2INC})
    set (XML_LIBS ${XML2LIB})
else ()
    find_package(LibXml2)
    if (LIBXML2_FOUND)
        set (XML_INCLUDES ${LIBXML2_INCLUDE_DIR})
        set (XML_LIBS ${LIBXML2_LIBRARIES})
    else ()
        set (XML_INCLUDES $ENV{XML2INC})
        set (XML_LIBS $ENV{XML2LIB})
    endif ()
endif ()
set (XML_INCLUDES ${XML_INCLUDES} CACHE PATH "XML2 header file path" FORCE)
set (XML_LIBS ${XML_LIBS} CACHE STRING "XML2 libraries" FORCE)

set (HDFEOS_GCTPINC ${HDFEOS_GCTPINC} CACHE PATH "HDFEOS GCTP header file path"
     FORCE)
set (HDFEOS_GCTPLIBS ${HDFEOS_GCTPLIB} CACHE STRING "HDFEOS GCTP libraries"
     FORCE)

set (ESPA_INCLUDES ${ESPA_INC} CACHE PATH "ESPA header file path" FORCE)
set (ESPA_LIBS ${ESPA_LIB} CACHE PATH "ESPA library path" FORCE)

enable_testing()

add_subdirectory(scripts)
add_subdirectory(provisional-st)
add_subdirectory(st_auxiliary_data)
