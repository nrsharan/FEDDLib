TRIBITS_ADD_EXECUTABLE(
	artery_dan_restart
    SOURCES main.cpp
)

# Restart test, two runs of the same case: the uninterrupted run
# (simulationParameters.xml) writes checkpoints at t = 0.05 and 0.1 to
# restartCheckpoints/; the restart (simulationParameters_restart.xml) continues
# from the checkpoint at t = 0.05, during the load ramp, and must reproduce the
# uninterrupted solution at t = 0.1 (restart during fiber reorientation is
# tested by sci_Cube_Restart).
# The case needs the AceGen elements (Interface2)
if(TPL_ENABLE_AceGENInterface)
  TRIBITS_ADD_ADVANCED_TEST(
      artery_dan_restart
      TEST_0 EXEC artery_dan_restart
          ARGS --simulationsParameters=simulationParameters.xml
          NUM_MPI_PROCS 4
      TEST_1 EXEC artery_dan_restart
          ARGS --simulationsParameters=simulationParameters_restart.xml
          NUM_MPI_PROCS 4
      OVERALL_NUM_MPI_PROCS 4
      TIMEOUT 3600
      FAIL_FAST
  )
endif()

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_artery_dan_restart
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES materialParameters_dan.xml preconditionerParameters_Chemistry.xml
	preconditionerParameters_Structure.xml simulationParameters.xml simulationParameters_restart.xml
	solverParameters.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES materialParameters_dan.xml preconditionerParameters_Chemistry.xml
	preconditionerParameters_Structure.xml simulationParameters.xml simulationParameters_restart.xml
	solverParameters.xml
	EXEDEPS artery_dan_restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_artery_dan_restart
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES Artery_dan_SCI.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES Artery_dan_SCI.mesh
	EXEDEPS artery_dan_restart
)
