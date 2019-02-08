add_test(NAME TestGridPoints COMMAND python3 -m unittest st-tests.TestGridPoints WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestAuxMerra COMMAND python3 -m unittest st-tests.TestAuxMerra WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestAuxGeos5 COMMAND python3 -m unittest st-tests.TestAuxGeos5 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestAuxNarr COMMAND python3 -m unittest st-tests.TestAuxNarr WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestModtranInput COMMAND python3 -m unittest st-tests.TestModtranInput WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestEmissivity COMMAND python3 -m unittest st-tests.TestEmissivity WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestEmissivityStdev COMMAND python3 -m unittest st-tests.TestEmissivityStdev WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestModtran COMMAND python3 -m unittest st-tests.TestModtran WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestAtmosParam COMMAND python3 -m unittest st-tests.TestAtmosParam WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestBuildSTData COMMAND python3 -m unittest st-tests.TestBuildSTData WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestDistToCloud COMMAND python3 -m unittest st-tests.TestDistToCloud WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestQA COMMAND python3 -m unittest st-tests.TestQA WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestConvert COMMAND python3 -m unittest st-tests.TestConvert WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME TestGenerateProducts COMMAND python3 -m unittest st-tests.TestGenerateProducts WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TEST TestAtmosParam
    PROPERTY ENVIRONMENT
    PATH=${CMAKE_CURRENT_BINARY_DIR}/../../src:$ENV{PATH})
set_tests_properties(TestGridPoints TestAuxMerra TestAuxGeos5 TestAuxNarr 
    TestModtranInput TestEmissivity TestEmissivityStdev TestModtran 
    TestAtmosParam TestBuildSTData TestDistToCloud TestQA TestConvert 
    TestGenerateProducts
    PROPERTIES RESOURCE_LOCK working_directory)
