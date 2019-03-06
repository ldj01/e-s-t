# Add this one first since it takes the longest
add_test(NAME GenerateProducts
         COMMAND python3 -m unittest st-tests.TestGenerateProducts)
add_test(NAME GridPoints COMMAND python3 -m unittest st-tests.TestGridPoints)
add_test(NAME AuxMerra COMMAND python3 -m unittest st-tests.TestAuxMerra)
add_test(NAME AuxGeos5 COMMAND python3 -m unittest st-tests.TestAuxGeos5)
add_test(NAME AuxNarr COMMAND python3 -m unittest st-tests.TestAuxNarr)
add_test(NAME ModtranInput
         COMMAND python3 -m unittest st-tests.TestModtranInput)
add_test(NAME Emissivity COMMAND python3 -m unittest st-tests.TestEmissivity)
add_test(NAME EmissivityStdev
         COMMAND python3 -m unittest st-tests.TestEmissivityStdev)
add_test(NAME Modtran COMMAND python3 -m unittest st-tests.TestModtran)
add_test(NAME AtmosParam COMMAND python3 -m unittest st-tests.TestAtmosParam)
add_test(NAME BuildSTData COMMAND python3 -m unittest st-tests.TestBuildSTData)
add_test(NAME DistToCloud COMMAND python3 -m unittest st-tests.TestDistToCloud)
add_test(NAME QA COMMAND python3 -m unittest st-tests.TestQA)
add_test(NAME Convert COMMAND python3 -m unittest st-tests.TestConvert)

set_property(TEST GridPoints AuxMerra AuxGeos5 AuxNarr ModtranInput 
             Emissivity EmissivityStdev Modtran AtmosParam BuildSTData 
             DistToCloud QA Convert GenerateProducts
             PROPERTY ENVIRONMENT
             PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/..:${CMAKE_CURRENT_SOURCE_DIR}:$ENV{PYTHONPATH}
             PATH=${CMAKE_BINARY_DIR}/provisional-st/src:${CMAKE_CURRENT_SOURCE_DIR}/..:$ENV{PATH})
