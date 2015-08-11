
cmake_minimum_required ( VERSION 2.8.12 )

install ( FILES __init__.py
                settings.py
          DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/config )
