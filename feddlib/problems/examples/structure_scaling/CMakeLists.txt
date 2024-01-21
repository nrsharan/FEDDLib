TRIBITS_ADD_EXECUTABLE(
	structure_scaling
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	structure_scaling
	NAME structure_scaling
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_structure_scaling
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml parametersProblemStructure.xml parametersProblemStructure_NH.xml parametersProblemStructure_SMC.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml parametersProblemStructure.xml parametersProblemStructure_NH.xml parametersProblemStructure_SMC.xml
	EXEDEPS structure_scaling
)

#TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_structure_scaling
#    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
#    SOURCE_FILES testFoam2.mesh
#    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
#    DEST_FILES testFoam2.mesh
#    EXEDEPS unsteadyNonLinElasticity
#)



