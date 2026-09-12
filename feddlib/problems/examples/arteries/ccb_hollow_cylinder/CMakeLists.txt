TRIBITS_ADD_EXECUTABLE(
	ccb_hollow_cylinder
    SOURCES main.cpp
)

# The case needs the AceGen elements (Interface2)
if(TPL_ENABLE_AceGENInterface)
  TRIBITS_ADD_TEST(
      ccb_hollow_cylinder
      NAME ccb_hollow_cylinder
      ARGS --simulationsParameters=simulationParameters_test.xml
      NUM_MPI_PROCS 4
  )
endif()

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_ccb_hollow_cylinder
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES materialParameters.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml simulationParameters_test.xml solverParameters.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES materialParameters.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml simulationParameters_test.xml solverParameters.xml
	EXEDEPS ccb_hollow_cylinder
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_ccb_hollow_cylinder
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/ccb_comparison
    SOURCE_FILES hollow_cylinder_p1.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES hollow_cylinder_p1.mesh
	EXEDEPS ccb_hollow_cylinder
)
