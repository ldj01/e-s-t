
cmake_minimum_required ( VERSION 2.8.12 )

find_package ( PythonInterp 2.7.0 REQUIRED )
find_package ( PythonLibs 2.7.0 REQUIRED )

install ( PROGRAMS do_lst.py
                   download_narr_aux_data.py
                   lst_extract_modtran_results.py
                   generate_aster_ndvi_and_b6_emis_mosaics.py
          DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )

