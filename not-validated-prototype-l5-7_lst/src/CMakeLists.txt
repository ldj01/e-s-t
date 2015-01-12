
cmake_minimum_required ( VERSION 2.8.12 )

include ( FindESPALibCommon.cmake )

find_package ( LibXml2 2.9.1 REQUIRED )
find_package ( ZLIB 1.2.8 REQUIRED )
find_package ( LibLZMA 5.1.2 REQUIRED )

find_library ( Math_Library m ) # We need the standard math library

include_directories ( ${LibESPA_INCLUDES} ${LIBXML2_INCLUDE_DIR} )

add_executable ( scene_based_lst input.c
                                 output.c
                                 utilities.c
                                 2d_array.c
                                 date.c
                                 misc.c
                                 split_filename.c
                                 build_modtran_input.c
                                 second_narr.c
                                 third_pixels_post.c
                                 scene_based_lst.c )

target_link_libraries ( scene_based_lst ${LibESPA_LIBRARIES}
                                        ${LIBXML2_LIBRARIES}
                                        ${ZLIB_LIBRARIES}
                                        ${LIBLZMA_LIBRARIES}
                                        ${Math_Library} )

install ( PROGRAMS scene_based_lst
          DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )

