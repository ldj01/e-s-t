
cmake_minimum_required ( VERSION 2.8.12 )

include ( FindESPALibCommon.cmake )

find_package ( LibXml2 2.9.1 REQUIRED )
find_package ( ZLIB 1.2.8 REQUIRED )
find_package ( LibLZMA 5.1.2 REQUIRED )

find_library ( Math_Library m ) # We need the standard math library

include_directories ( ${LibESPA_INCLUDES} ${LIBXML2_INCLUDE_DIR} )

add_executable ( lst lst.c
                     input.c
                     output.c
                     utilities.c
                     2d_array.c
                     date.c
                     get_args.c
                     split_filename.c
                     build_points.c
                     build_modtran_input.c
                     calculate_point_atmospheric_parameters.c
                     calculate_pixel_atmospheric_parameters.c )

target_link_libraries ( lst ${LibESPA_LIBRARIES}
                            ${LIBXML2_LIBRARIES}
                            ${ZLIB_LIBRARIES}
                            ${LIBLZMA_LIBRARIES}
                            ${Math_Library} )

install ( TARGETS lst
          DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )

