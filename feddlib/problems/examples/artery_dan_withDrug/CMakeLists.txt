TRIBITS_ADD_EXECUTABLE(
	 artery_dan_withDrug
    SOURCES main.cpp
)

TRIBITS_ADD_TEST(
     artery_dan_withDrug
    NAME  artery_dan_withDrug
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_ artery_dan_withDrug
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES materialParameters_dan.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml solverParameters.xml slurm-script.sh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES materialParameters_dan.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml solverParameters.xml slurm-script.sh
	EXEDEPS  artery_dan_withDrug
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_ artery_dan_withDrug
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES  artery_dan_SCI.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES  artery_dan_SCI.mesh
	EXEDEPS  artery_dan_withDrug
)
