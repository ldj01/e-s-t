
cmake_minimum_required ( VERSION 2.8.12 )

install ( FILES modtran_head.txt
                modtran_tail.txt
                narr_coordinates.txt
                std_mid_lat_summer_atmos.txt
          DESTINATION ${CMAKE_INSTALL_PREFIX}/static_data )

