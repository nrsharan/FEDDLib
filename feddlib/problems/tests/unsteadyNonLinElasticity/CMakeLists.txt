TRIBITS_ADD_EXECUTABLE(
	unsteadyNonLinElasticityAssFE
	unsteadyNonLinElasticityAssFE
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	unsteadyNonLinElasticityAssFE
	NAME unsteadyNonLinElasticityAssFE
	unsteadyNonLinElasticityAssFE
	NAME unsteadyNonLinElasticityAssFE
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadyNonLinElasticityAssFE
TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadyNonLinElasticityAssFE
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml parametersProblemStructure.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml parametersProblemStructure.xml
	EXEDEPS unsteadyNonLinElasticityAssFE
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyNonLinElasticityAssFE
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES artery_h_0_8_L_0_25.mesh cube_h_1.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES artery_h_0_8_L_0_25.mesh cube_h_1.mesh
    EXEDEPS unsteadyNonLinElasticityAssFE
)
