set (SCRIPTS st_determine_grid_points.py st_build_modtran_input.py
             st_extract_auxiliary_narr_data.py
             st_generate_distance_to_cloud.py st_generate_products.py
             st_generate_qa.py st_run_modtran.py
             estimate_landsat_emissivity.py
             estimate_landsat_emissivity_stdev.py st_convert_bands.py)
set (SCRIPT_IMPORTS build_st_data.py st_exceptions.py st_grid_points.py
                    st_utilities.py emissivity_utilities.py)

install(PROGRAMS ${SCRIPTS} DESTINATION bin)
install(FILES ${SCRIPT_IMPORTS} DESTINATION bin)

add_subdirectory(unit-tests)
