# Set version info using Git tags (if available) or a stored version.txt file.
# Since versioning info isn't always in the standard "x.x.x" format, we need
# to build up the version based on what is available by delimiting on
# periods as we identify major, minor, and patch components.
execute_process(COMMAND ${PROJECT_SOURCE_DIR}/get_version.pl
                        ${PROJECT_SOURCE_DIR}
                OUTPUT_VARIABLE VERSION_STRING
                OUTPUT_STRIP_TRAILING_WHITESPACE)

string(REPLACE "." ";" VERSION_LIST ${VERSION_STRING})
list(LENGTH VERSION_LIST VERSION_LENGTH)
list(GET VERSION_LIST 0 VERSION_MAJOR)
set (VERSION ${VERSION_MAJOR})
set (SHORT_VERSION "${VERSION_MAJOR}")
set (VERSION_MINOR 0)
set (VERSION_PATCH 0)
set (VERSION_SHA1 0)
if (${VERSION_LENGTH} GREATER 1)
    list(GET VERSION_LIST 1 VERSION_MINOR)
    string(APPEND VERSION ".${VERSION_MINOR}")
    string(APPEND SHORT_VERSION ".${VERSION_MINOR}")
    if (${VERSION_LENGTH} GREATER 2)
        list(GET VERSION_LIST 2 VERSION_PATCH)
        string(APPEND VERSION ".${VERSION_PATCH}")
        string(APPEND SHORT_VERSION ".${VERSION_PATCH}")
        if (${VERSION_LENGTH} GREATER 3)
            list(GET VERSION_LIST 3 VERSION_SHA1)
            string(APPEND VERSION ".${VERSION_SHA1}")
        endif ()
    endif ()
endif ()

message(STATUS "Package version: ${VERSION}")

set (${PROJECT_NAME}_VERSION_MAJOR ${VERSION_MAJOR})
set (${PROJECT_NAME}_VERSION_MINOR ${VERSION_MINOR})
set (${PROJECT_NAME}_VERSION_PATCH ${VERSION_PATCH})

set (PACKAGE ${PROJECT_NAME})

set(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
set(CPACK_PACKAGE_VERSION ${SHORT_VERSION})
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${VERSION}")
set(CPACK_SOURCE_IGNORE_FILES "/build/;/.git/;/.gz/;~$;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)
