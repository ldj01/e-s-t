file(GLOB STATIC_DATA "*.txt")

install(FILES ${STATIC_DATA} DESTINATION static_data)
