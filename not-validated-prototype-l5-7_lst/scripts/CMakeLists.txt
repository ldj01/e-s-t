
cmake_minimum_required ( VERSION 2.8.12 )

find_package ( PythonInterp 2.7.0 REQUIRED )
find_package ( PythonLibs 2.7.0 REQUIRED )

add_subdirectory ( config )

install ( PROGRAMS l5_7_lst.py
                   l5_7_build_lst_data.py
                   l5_7_estimate_landsat_emissivity.py
                   download_narr_aux_data.py
                   lst_extract_modtran_results.py
                   lst_prototype.py
                   update_narr_aux_data_from_CISL_RDA.py
          DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )

install ( FILES lst_utilities.py
          DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )
