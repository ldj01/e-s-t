
cmake_minimum_required ( VERSION 2.6.0 )

include ( FindESPALibCommon.cmake )

find_package ( LibXml2 2.9.1 REQUIRED )
find_package ( ZLIB 1.2.8 REQUIRED )
find_package ( LibLZMA 5.1.2 REQUIRED )

find_library ( Math_Library m ) # We need the standard math library

include_directories ( ${LibESPA_INCLUDES} ${LIBXML2_INCLUDE_DIR} )

add_executable ( scene_based_lst error.c input.c 2d_array.c scene_based_lst.c
                                 first_files.c second_narr.c
                                 third_pixels_post.c date.c split_filename.c
                                 misc.c output.c )

target_link_libraries ( scene_based_lst ${LibESPA_LIBRARIES}
                                        ${LIBXML2_LIBRARIES}
                                        ${ZLIB_LIBRARIES}
                                        ${LIBLZMA_LIBRARIES}
                                        ${Math_Library} )

install ( PROGRAMS scene_based_lst
          DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )

