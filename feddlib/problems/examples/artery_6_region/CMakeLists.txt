TRIBITS_ADD_EXECUTABLE(
	artery_6_region
    SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    artery_6_region
    NAME artery_6_region
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_artery_6_region
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES materialParameters.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml 
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES materialParameters.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml
	EXEDEPS artery_6_region
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_artery_6_region
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES artery_h_025_L_025.mesh artery_h_0_8_L_0_25.mesh plaque_short_6_k.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES artery_h_025_L_025.mesh artery_h_0_8_L_0_25.mesh plaque_short_6_k.mesh
	EXEDEPS artery_6_region
)
