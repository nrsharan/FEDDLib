SET(Mesh_HEADERS
    Mesh/gmsh_io.hpp
	Mesh/MeshUnstructured.hpp
	Mesh/MeshUnstructured_decl.hpp
	Mesh/MeshUnstructured_def.hpp
	Mesh/MeshFactory.hpp
    Mesh/MeshFactory_decl.hpp
	Mesh/MeshFactory_def.hpp
	Mesh/MeshFileReader.hpp
	Mesh/MeshInterface.hpp
    Mesh/MeshInterface_decl.hpp
    Mesh/MeshInterface_def.hpp
    Mesh/MeshPartitioner.hpp
    Mesh/MeshPartitioner_decl.hpp
    Mesh/MeshPartitioner_def.hpp
    Mesh/AABBTree.hpp
    Mesh/AABBTree_decl.hpp
    Mesh/AABBTree_def.hpp
CACHE INTERNAL "")

SET(Mesh_SOURCES
    Mesh/gmsh_io.cpp
    Mesh/MeshUnstructured.cpp
    Mesh/MeshFactory.cpp
    Mesh/MeshPartitioner.cpp
    Mesh/MeshFileReader.cpp
    Mesh/MeshInterface.cpp
    Mesh/AABBTree.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${Mesh_HEADERS})
  STRING(REPLACE "Mesh/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/core/Mesh"
  COMPONENT ${PACKAGE_NAME}
)
