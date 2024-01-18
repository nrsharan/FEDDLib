TRIBITS_ADD_EXECUTABLE(
	sci_scaling
    SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    sci_scaling
    NAME sci_scaling
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_sciscaling
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecChem.xml parametersPrec.xml parametersPrecTeko.xml parametersPrecBlock.xml parametersProblemSCI.xml parametersProblemStructure.xml parametersSolverSCI.xml parametersPrecStructure.xml parametersProblemStructure_plaque_full.xml sbatchStrongScaling100Mono.sh parametersPrecCE.xml
	#sbatchWeakScaling216Mono.sh sbatchWeakScaling1000Mono.sh sbatchWeakScaling4096Mono.sh sbatchStrongScaling100Mono.sh sbatchStrongScaling200Mono.sh sbatchStrongScaling400Mono.sh sbatchStrongScaling800Mono.sh 
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecChem.xml parametersPrec.xml parametersPrecTeko.xml parametersPrecBlock.xml parametersProblemSCI.xml parametersProblemStructure.xml parametersSolverSCI.xml parametersPrecStructure.xml parametersProblemStructure_plaque_full.xml sbatchStrongScaling100Mono.sh parametersPrecCE.xml
	# sbatchWeakScaling216Mono.sh sbatchWeakScaling1000Mono.sh sbatchWeakScaling4096.sh sbatchStrongScaling100Mono.sh sbatchStrongScaling200Mono.sh sbatchStrongScaling400Mono.sh sbatchStrongScaling800Mono.sh 
	EXEDEPS sci_scaling
)
# cube_h_0_05.mesh cube_h_0_1.mesh artery_h_025.mesh artery_h_0_1.mesh structure_realistic_artery.mesh
TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_sciscaling
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES  cube_h_1.mesh artery_h_025_L_025.mesh artery_h_0_8_L_0_25.mesh artery_h_015_L_075.mesh Artery_Small.mesh Artery_Long.mesh plaque_4k.mesh arterie_with_plaque_initma_media_adventitia.mesh arterie_with_plaque_thin_initma_media_adventitia.mesh plaque_short_6_k.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES cube_h_1.mesh  artery_h_025_L_025.mesh artery_h_0_8_L_0_25.mesh artery_h_015_L_075.mesh Artery_Small.mesh Artery_Long.mesh plaque_4k.mesh arterie_with_plaque_initma_media_adventitia.mesh arterie_with_plaque_thin_initma_media_adventitia.mesh plaque_short_6_k.mesh
	EXEDEPS sci_scaling 
)
