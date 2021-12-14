TRIBITS_ADD_EXECUTABLE(
	steadyNonlinearReactionDiffusion
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	steadyNonlinearReactionDiffusion
	NAME steadyNonlinearReactionDiffusion
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_steadyNonlinearReactionDiffusion
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml 
	EXEDEPS steadyNonlinearReactionDiffusion
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_steadyNonlinearReactionDiffusion
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES square.mesh
	EXEDEPS steadyNonlinearReactionDiffusion
)

