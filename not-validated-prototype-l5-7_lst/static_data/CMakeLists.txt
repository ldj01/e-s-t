
cmake_minimum_required ( VERSION 2.8.12 )

install ( FILES modtran_head.txt
                modtran_tail.txt
                narr_coordinates.txt
                std_mid_lat_summer_atmos.txt
                L5_Brightness_Temperature_LUT.txt
                L5_Spectral_Response.txt
                L7_Brightness_Temperature_LUT.txt
                L7_Spectral_Response.txt
                L8_Spectral_Response.txt
          DESTINATION ${CMAKE_INSTALL_PREFIX}/static_data )

