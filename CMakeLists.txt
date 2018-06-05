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
        set (XML_LIBS -L$ENV{XML2LIB} -lxml2)
    endif ()
endif ()
set (XML_INCLUDES ${XML_INCLUDES} CACHE PATH "XML2 header file path" FORCE)
set (XML_LIBS ${XML_LIBS} CACHE STRING "XML2 libraries" FORCE)

if (NOT GCTP_INCLUDES)
    set (GCTP_INCLUDES $ENV{HDFEOS_GCTPINC})
endif ()
if (NOT GCTP_LIBS)
    set (GCTP_LIBS "-L$ENV{HDFEOS_GCTPLIB} -lGctp -lm")
endif ()
set (GCTP_INCLUDES ${GCTP_INCLUDES} CACHE PATH "GCTP header file path" FORCE)
set (GCTP_LIBS ${GCTP_LIBS} CACHE STRING "GCTP libraries" FORCE)

set (ESPA_INCLUDES $ENV{ESPAINC} CACHE PATH "ESPA header file path" FORCE)
set (ESPA_LIBS $ENV{ESPALIB} CACHE PATH "ESPA library path" FORCE)

enable_testing()

add_custom_target(rit-aux
                  COMMAND "${CMAKE_COMMAND}"
                          -DCMAKE_INSTALL_COMPONENT=rit-aux
                          -P "${CMAKE_BINARY_DIR}/cmake_install.cmake")

add_subdirectory(scripts)
add_subdirectory(provisional-st)
add_subdirectory(st_auxiliary_data)
