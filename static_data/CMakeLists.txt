
cmake_minimum_required ( VERSION 2.6.0 )

find_package ( PythonInterp 2.7.0 REQUIRED )
find_package ( PythonLibs 2.7.0 REQUIRED )

install ( FILES modtran_head.txt modtran_tail.txt narr_coordinates.txt
          DESTINATION ${CMAKE_INSTALL_PREFIX}/static_data )

