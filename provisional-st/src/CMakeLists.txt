configure_file(const.h.in const.h)

include_directories(${ESPA_INCLUDES})
include_directories(${XML_INCLUDES})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

link_directories(${ESPA_LIBS})

file(GLOB SOURCES "*.c")

add_executable(st_atmospheric_parameters ${SOURCES})
target_link_libraries(st_atmospheric_parameters -lespa_raw_binary -lespa_common
                      -lespa_format_conversion
                      ${GCTP2_LIBS} ${XML_LIBS} m)
install(TARGETS st_atmospheric_parameters DESTINATION bin)
