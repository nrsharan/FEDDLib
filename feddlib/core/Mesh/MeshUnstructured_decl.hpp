#ifndef MESHUNSTRUCTURED_decl_hpp
#define MESHUNSTRUCTURED_decl_hpp

#include "feddlib/core/Utils/FEDDUtils.hpp"
#include "MeshInterface.hpp"
#include "MeshFileReader.hpp"
#include "feddlib/core/FE/EdgeElements.hpp"
#include "feddlib/core/FE/TriangleElements.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/FE/Elements.hpp"
#include "feddlib/core/Mesh/AABBTree.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
/*!
 Declaration of MeshUnstructured
 
 @brief  MeshUnstructured
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
    
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class MeshUnstructured{
    
public:
    typedef MeshUnstructured<SC,LO,GO,NO> Mesh_Type;
    typedef Teuchos::RCP<MeshUnstructured<SC,LO,GO,NO> > MeshPtr_Type;

    typedef std::vector<MeshPtr_Type> MeshPtrArray_Type;

	typedef Elements Elements_Type;
    typedef Teuchos::RCP<Elements_Type> ElementsPtr_Type;
    
    typedef Teuchos::RCP<Teuchos::Comm<int> > CommPtr_Type;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > CommConstPtr_Type;
    typedef const CommConstPtr_Type CommConstPtrConst_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef Teuchos::RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtrConst_Type;


    typedef EdgeElements EdgeElements_Type;
    typedef Teuchos::RCP<EdgeElements_Type> EdgeElementsPtr_Type;

    typedef SurfaceElements SurfaceElements_Type;
    typedef Teuchos::RCP<SurfaceElements_Type> SurfaceElementsPtr_Type;
    
    typedef MeshInterface<SC,LO,GO,NO> MeshInterface_Type;
    typedef Teuchos::RCP<MeshInterface_Type> MeshInterfacePtr_Type;
   
    typedef Teuchos::OrdinalTraits<LO> OTLO;

	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
	typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
	typedef MultiVector<LO,LO,GO,NO> MultiVectorLO_Type;
	typedef Teuchos::RCP<MultiVectorLO_Type> MultiVectorLOPtr_Type;
	typedef Teuchos::RCP<const MultiVector_Type> MultiVectorPtrConst_Type;

    typedef AABBTree<SC,LO,GO,NO> AABBTree_Type;
    typedef Teuchos::RCP<AABBTree_Type > AABBTreePtr_Type;

	typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;
    
    MeshUnstructured();
    
    MeshUnstructured( CommConstPtr_Type comm, int volumeID=10 );
    
    //~MeshUnstructured();

    virtual ~MeshUnstructured();
    
//    virtual vec2D_int_ptr_Type getElements();

   // virtual void dummy() = 0;

    void setParameterList( ParameterListPtr_Type& pL );

    ParameterListConstPtr_Type getParameterList( ) const;
    
    vec_int_ptr_Type getElementsFlag() const;
    
    MapConstPtr_Type getMapUnique() const;
    
    MapConstPtr_Type getMapRepeated() const;

    MapConstPtr_Type getMapUniqueP2() const;
    
    MapConstPtr_Type getMapRepeatedP2() const;
    
    MapConstPtr_Type getElementMap();
	
    MapConstPtr_Type getEdgeMap(); // Edge Map
    
    vec2D_dbl_ptr_Type getPointsRepeated() const;

    vec2D_dbl_ptr_Type getPointsUnique() const;
    
    vec_int_ptr_Type getBCFlagRepeated() const;
    
    vec_int_ptr_Type getBCFlagUnique() const;
    
    vec2D_int_ptr_Type getElements() {
        vec2D_int_ptr_Type tmp;
        return tmp;
    };
    
    ElementsPtr_Type getElementsC();

    ElementsPtr_Type getSurfaceElements();
    
    int getDimension();
    
    GO getNumElementsGlobal();
    
    LO getNumElements();
    
    LO getNumPoints(std::string type="Unique");
    
    int getOrderElement();

    CommConstPtr_Type getComm(){return comm_;};

	// NEW
	int setStructuredMeshFlags(int flags){return 0;};
    
    void setElementFlags(std::string type="");
    
    void setReferenceConfiguration();

 	tuple_intint_Type getRankRange() const {return rankRange_;};
    
    void deleteSurfaceElements(){ surfaceElements_.reset(); };
    
    void setMeshFileName(string meshFileName, string delimiter);

 	//int determineFlagP2( FiniteElement& fe, LO p1ID, LO p2ID, vec2D_int_Type& permutation );
    
    int determineFlagP2( LO p1ID, LO p2ID,  LO localEdgeID, vec2D_LO_Type& markedPoint );
	// ----------------------------------------------------------------
     
    void getTriangles(int vertex1ID, int vertex2ID, vec_int_Type &vertices3ID);

    SurfaceElementsPtr_Type getSurfaceTriangleElements(){return surfaceTriangleElements_;};
    
    void findSurfaces( const vec_int_Type& elementNodeList, vec_int_Type numbering,  vec2D_int_Type& localSurfaceNodeList_vec, vec_int_Type& locSurfaces, bool critical = false );

    void findEdges( const vec_int_Type& elementNodeList, vec_int_Type numbering,  vec2D_int_Type& localEdgeNodeList_vec, vec_int_Type& locEdges);
    
    MeshInterfacePtr_Type getMeshInterface();        
    
    void setEdgeElements( EdgeElementsPtr_Type edgeElements ){ edgeElements_ = edgeElements; };
    
    EdgeElementsPtr_Type getEdgeElements( ){ return edgeElements_; };
    
    ElementsPtr_Type getSurfaceEdgeElements(){return surfaceEdgeElements_;};
    
    int getSurfaceElementOrder(){return surfaceElementOrder_;};
    
    int getEdgeElementOrder(){return edgesElementOrder_;};
    
    int getNumGlobalNodes(){return numNodes_;};
    
    // Creates an AABBTree from own vertice- and elementlist.
    void create_AABBTree();
    
    vec_int_ptr_Type  findElemsForPoints(vec2D_dbl_ptr_Type query_points);
    
    vec_dbl_Type getBaryCoords(vec_dbl_Type point, int element);
    
    bool isPointInElem(vec_dbl_Type point, int element);

	// additional functions to complete mesh class.
	void buildEdgeMap();

	void updateElementsOfEdgesLocalAndGlobal();

	void assignEdgeFlags();
  
    /* ###################################################################### */
    
    MeshInterfacePtr_Type meshInterface_;
        
    int volumeID_;
    /* ###################################################################### */
 	string meshFileName_;
    string delimiter_;

    int elementOrder_;
    int surfaceElementOrder_;
    int edgesElementOrder_;
    int numElements_;
    int numSurfaces_;
    int numEdges_;
    int numNodes_;

    int                     dim_;
    long long               numElementsGlob_;

    std::string 			FEType_;
    MapPtr_Type             mapUnique_;
    MapPtr_Type 			mapRepeated_;
    vec2D_dbl_ptr_Type		pointsRep_;
    vec2D_dbl_ptr_Type 		pointsUni_;
    vec_int_ptr_Type 		bcFlagRep_;
    vec_int_ptr_Type		bcFlagUni_;

    ElementsPtr_Type        surfaceElements_;

    ElementsPtr_Type        elementsC_;

	// NEW:
 	EdgeElementsPtr_Type edgeElements_;    
    ElementsPtr_Type surfaceEdgeElements_;
	SurfaceElementsPtr_Type surfaceTriangleElements_;


    MapPtr_Type				elementMap_;
    MapPtr_Type				edgeMap_;


    CommConstPtr_Type  comm_;
    
    vec2D_dbl_ptr_Type		pointsRepRef_; // Repeated Referenzkonfiguration
    vec2D_dbl_ptr_Type		pointsUniRef_; // Unique Referenzkonfiguration
    
    //vec_int_ptr_Type 		elementFlag_;

    MapPtr_Type mapUniqueP2Map_;
    MapPtr_Type mapRepeatedP2Map_;
    
    ParameterListPtr_Type pList_;

    AABBTreePtr_Type AABBTree_;
    
    tuple_intint_Type rankRange_;

private:
    
     
};
}
#endif
