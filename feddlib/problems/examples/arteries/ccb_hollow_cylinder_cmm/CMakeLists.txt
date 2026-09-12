TRIBITS_ADD_EXECUTABLE(
	ccb_hollow_cylinder_cmm
    SOURCES main.cpp
)

# The case needs the AceGen elements (Interface2)
if(TPL_ENABLE_AceGENInterface)
  TRIBITS_ADD_TEST(
      ccb_hollow_cylinder_cmm
      NAME ccb_hollow_cylinder_cmm
      ARGS --simulationsParameters=simulationParameters_test.xml
      NUM_MPI_PROCS 4
  )
endif()

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_ccb_hollow_cylinder_cmm
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES materialParameters.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml simulationParameters_test.xml solverParameters.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES materialParameters.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml simulationParameters_test.xml solverParameters.xml
	EXEDEPS ccb_hollow_cylinder_cmm
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_ccb_hollow_cylinder_cmm
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/ccb_comparison
    SOURCE_FILES hollow_cylinder_p1.mesh hollow_cylinder_p2_override.txt
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES hollow_cylinder_p1.mesh hollow_cylinder_p2_override.txt
	EXEDEPS ccb_hollow_cylinder_cmm
)
