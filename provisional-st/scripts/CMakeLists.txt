file(GLOB SCRIPTS "*.py")

install(PROGRAMS ${SCRIPTS} DESTINATION bin)

add_subdirectory(unit-tests)
