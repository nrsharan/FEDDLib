TRIBITS_ADD_EXECUTABLE(
	artery_narula_withDrug
    SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    artery_narula_withDrug
    NAME artery_narula_withDrug
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_artery_narula_withDrug
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES materialParameters_narula.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml solverParameters.xml slurm-script.sh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES materialParameters_narula.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml solverParameters.xml slurm-script.sh
	EXEDEPS artery_narula_withDrug
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_artery_narula_withDrug
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES Artery_narula_etal-coronoary_tcfa_medium_SCI.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES Artery_narula_etal-coronoary_tcfa_medium_SCI.mesh
	EXEDEPS artery_narula_withDrug
)
