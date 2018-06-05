file(GLOB STATIC_DATA "*.txt")
list(FILTER STATIC_DATA EXCLUDE REGEX "CMakeLists.txt")

install(FILES ${STATIC_DATA} DESTINATION static_data)
