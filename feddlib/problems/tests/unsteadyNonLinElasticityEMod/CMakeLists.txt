TRIBITS_ADD_EXECUTABLE(
	unsteadyNonLinElasticityEMod
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	unsteadyNonLinElasticityEMod
	NAME unsteadyNonLinElasticityEMod
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadyNonLinElasticityEMod
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS unsteadyNonLinElasticity
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyNonLinElasticityEMod
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES testFoam2.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES testFoam2.mesh
    EXEDEPS unsteadyNonLinElasticity
)
