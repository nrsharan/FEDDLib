TRIBITS_ADD_EXECUTABLE(
	test_flowrate
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_test_flowrate
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersProblem.xml plistProblemLaplace.xml plistPrecLaplace.xml plistSolverLaplace.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersProblem.xml plistProblemLaplace.xml plistPrecLaplace.xml plistSolverLaplace.xml
	EXEDEPS test_flowrate
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_test_flowrate
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES fluid_safe.mesh 
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES fluid_safe.mesh 
	EXEDEPS test_flowrate
)
