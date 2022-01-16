#ifndef MESHFACTORY_decl_hpp
#define MESHFACTORY_decl_hpp

#include "feddlib/core/Utils/FEDDUtils.hpp"
#include "MeshUnstructured.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "MeshInterface.hpp"
#include "MeshFileReader.hpp"

/*!
Defintion of Mesh

@brief  MeshFactory
@author Christian Hochmuth
@version 1.0
@copyright CH
*/

using namespace std;

namespace FEDD {
    
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class MeshFactory : public MeshUnstructured<SC,LO,GO,NO> {
    
public:

    typedef MeshUnstructured<SC,LO,GO,NO> Mesh_Type;
    typedef Teuchos::RCP<MeshUnstructured<SC,LO,GO,NO> > MeshPtr_Type;
    typedef std::vector<MeshPtr_Type> MeshUnstrPtrArray_Type;

 	typedef MeshFactory<SC,LO,GO,NO> MeshFactory_Type;
    typedef Teuchos::RCP<MeshFactory<SC,LO,GO,NO> > MeshFactoryPtr_Type;
    typedef std::vector<MeshFactoryPtr_Type> MeshFactoryPtrArray_Type;

    typedef typename Mesh_Type::CommPtr_Type CommPtr_Type;
    typedef typename Mesh_Type::CommConstPtr_Type CommConstPtr_Type;
    typedef typename Mesh_Type::Elements_Type Elements_Type;
    typedef typename Mesh_Type::ElementsPtr_Type ElementsPtr_Type;    
    
    typedef EdgeElements EdgeElements_Type;
    typedef Teuchos::RCP<EdgeElements_Type> EdgeElementsPtr_Type;

    typedef SurfaceElements SurfaceElements_Type;
    typedef Teuchos::RCP<SurfaceElements_Type> SurfaceElementsPtr_Type;
    
    typedef MeshInterface<SC,LO,GO,NO> MeshInterface_Type;
    typedef Teuchos::RCP<MeshInterface_Type> MeshInterfacePtr_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef typename Map_Type::MapPtr_Type MapPtr_Type;
    typedef typename Map_Type::MapConstPtr_Type MapConstPtr_Type;

    typedef Teuchos::OrdinalTraits<LO> OTLO;

	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
	typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
	typedef MultiVector<LO,LO,GO,NO> MultiVectorLO_Type;
	typedef Teuchos::RCP<MultiVectorLO_Type> MultiVectorLOPtr_Type;
	typedef Teuchos::RCP<const MultiVector_Type> MultiVectorPtrConst_Type;


    /* ###################################################################### */
    //
    MeshFactory();

    MeshFactory( CommConstPtr_Type comm,int volumeID=10 );
    
    ~MeshFactory();
    
    // virtual void dummy() {};
    /*!
     Delete all member variables
     */
    //void deleteData();
      
	// Building P2-Mesh
    
   	void buildP2ofP1MeshEdge( MeshPtr_Type meshP1 );

    void setP2SurfaceElements( MeshPtr_Type meshP1 );
    
    void setSurfaceP2( FiniteElement &feP2, const FiniteElement &surfFeP1, const vec2D_int_Type &surfacePermutation, int dim );
    
    vec_int_Type reorderP2SurfaceIndices( vec_int_Type& additionalP2IDs, vec_int_Type& index, bool track=false);
    
    void getLocalSurfaceIndices( vec2D_int_Type& surfacePermutation, int surfaceElementOrder );
    
    void getEdgeCombinations( vec2D_int_Type& edgeCombinations );
        
    void determinePositionInElementP2( vec_int_Type& positions, vec_GO_Type& elementsGlobalOfEdge, LO p1ID, LO p2ID, MeshPtr_Type meshP1 );   
	// ---------------
	// Reading Mesh
	void readMeshSize();
    
    void readMeshEntity(string entityType);
	// ---------------

	void buildMeshInterfaceParallelAndDistance( MeshPtr_Type mesh, vec_int_Type flag_vec, vec_dbl_ptr_Type &distancesToInterface );

    void partitionInterface();
	// ---------------
    void moveMesh( MultiVectorPtr_Type displacementUnique, MultiVectorPtr_Type displacementRepeated );

	// ---------------
	// building structured meshes
    void setGeometry2DRectangle(std::vector<double> coordinates, double l, double h);
    
    void setGeometry3DBox(std::vector<double> coordinates, double l, double w, double h);
    
    void buildMesh2D(std::string FEType,
                    int N,
                    int M,
                    int numProcsCoarseSolve=0,
                    std::string underlyingLib="Tpetra");

    void buildMesh2DTPM(std::string FEType,
                        int N,
                        int M,
                        int numProcsCoarseSolve=0,
                        std::string underlyingLib="Tpetra");
    
    void buildMesh3D(std::string FEType,
                     int N,
                     int M,
                     int numProcsCoarseSolve=0,
                     std::string underlyingLib="Tpetra");
    
    void buildMesh2DBFS(std::string FEType,
                     int N,
                     int M,
                     int numProcsCoarseSolve=0,
                     std::string underlyingLib="Tpetra");

    void buildMesh3DBFS(std::string FEType,
                        int N,
                        int M,
                        int numProcsCoarseSolve=0,
                        std::string underlyingLib="Tpetra");
    
    void buildP1_Disc_Q2_3DBFS(int N,
                            int M,
                            int numProcsCoarseSolve,
                            std::string underlyingLib);

    void buildP1_Disc_Q2_3DCube(int N,
                               int M,
                               int numProcsCoarseSolve,
                               std::string underlyingLib);

    void build3DQ1Cube(int N,
                       int M,
                       int numProcsCoarseSolve,
                       std::string underlyingLib );

    void build3DQ2Cube( int N,
                        int M,
                        int numProcsCoarseSolve,
                        std::string underlyingLib );

    void build3DQ2_20Cube(int N,
                          int M,
                          int numProcsCoarseSolve,
                          std::string underlyingLib ); 


    void build3DQ2BFS( int N,
                       int M,
                       int numProcsCoarseSolve,
                       std::string underlyingLib );
    
    void buildMesh2DMiniTPM(std::string FEType,
                            int N,
                            int M,
                            int numProcsCoarseSolve=0,                            
                            std::string underlyingLib="Tpetra" );

    void buildSurfaceLinesSquareMiniTPM( string feType );
    
    void setRankRange( int numProcsCoarseSolve );
    
    void buildElementsClass( vec2D_int_ptr_Type elements, vec_int_ptr_Type elementFlag = Teuchos::null );
    
    void buildSurfaceLinesSquare();
    
    GO globalID_Q2_20Cube(int r, int s , int t, int &rr, int off_x, int off_y, int off_z, int M, int N,
                          GO nmbPoints_oneDirFull, GO nmbPoints_oneDirMid);
    
    void setStructuredMeshFlags(int flagsOption,string FEType="P1");
    
    void buildElementMap();

	void buildEdgeElements();

	void assignEdgeFlags();

	// ----------------------
    MapPtr_Type mapUniqueP2Map_;
    MapPtr_Type mapRepeatedP2Map_;

    std::vector<double> coorRec;
    double 				length;
    double		 		height;
    double 				width;

    //AABBTreePtr_Type AABBTree_;
    
    /* ###################################################################### */
private:
    void readSurfaces();
    
    void readLines();
    
    void readElements();
    
    void readNodes();

};
}

#endif
