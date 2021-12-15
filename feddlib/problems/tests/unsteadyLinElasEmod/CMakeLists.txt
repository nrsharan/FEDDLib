TRIBITS_ADD_EXECUTABLE(
	unsteadyLinElasEmod
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	unsteadyLinElasEmod
	NAME unsteadyLinElasEmod
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 2
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadyLinElasEmod
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS unsteadyLinElasEmod
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyLinElasEmod
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES Rohrhalb.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES Rohrhalb.mesh
	EXEDEPS unsteadyLinElasEmod
)
