TRIBITS_ADD_EXECUTABLE(
	simpleTest
    SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    simpleTest
    NAME simpleTest
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_simpleTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES materialParameters.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml solverParameters.xml slurm-script.sh slurm-script_elysium.sh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES materialParameters.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml solverParameters.xml slurm-script.sh slurm-script_elysium.sh
	EXEDEPS simpleTest
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_simpleTest
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES fsi_solid_2mm.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES fsi_solid_2mm.mesh
	EXEDEPS simpleTest
)
