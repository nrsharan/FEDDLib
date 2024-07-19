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
	SOURCE_FILES materialParameters_kim_guzman.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml solverParameters.xml slurm-script.sh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES materialParameters_kim_guzman.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml solverParameters.xml slurm-script.sh
	EXEDEPS artery_6_region
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_artery_6_region
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES Artery_kim_guzman_SCI.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES Artery_kim_guzman_SCI.mesh
	EXEDEPS artery_6_region
)
