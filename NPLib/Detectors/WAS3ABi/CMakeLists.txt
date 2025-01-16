add_custom_command(OUTPUT TWAS3ABiPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TWAS3ABiPhysics.h TWAS3ABiPhysicsDict.cxx TWAS3ABiPhysics.rootmap libNPWAS3ABi.dylib DEPENDS TWAS3ABiPhysics.h)
add_custom_command(OUTPUT TWAS3ABiDataDict.cxx COMMAND ../../scripts/build_dict.sh TWAS3ABiData.h TWAS3ABiDataDict.cxx TWAS3ABiData.rootmap libNPWAS3ABi.dylib DEPENDS TWAS3ABiData.h)
add_library(NPWAS3ABi SHARED TWAS3ABiSpectra.cxx TWAS3ABiData.cxx TWAS3ABiPhysics.cxx TWAS3ABiDataDict.cxx TWAS3ABiPhysicsDict.cxx )
target_link_libraries(NPWAS3ABi ${ROOT_LIBRARIES} NPCore) 
install(FILES TWAS3ABiData.h TWAS3ABiPhysics.h TWAS3ABiSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

