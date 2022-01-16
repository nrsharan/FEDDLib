#ifndef MESHUNSTRUCTURED_def_hpp
#define MESHUNSTRUCTURED_def_hpp
#include "MeshUnstructured_decl.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector_def.hpp"

/*!
 Definition of MeshUnstructured
 
 @brief  MeshUnstructured
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_MIN;
using Teuchos::outArg;

namespace FEDD {

template <class SC, class LO, class GO, class NO>
MeshUnstructured<SC,LO,GO,NO>::MeshUnstructured():
meshInterface_(),
numElementsGlob_(0),
mapUnique_(),
mapRepeated_(),
pointsRep_(),
pointsUni_(),
bcFlagRep_(),
bcFlagUni_(),
surfaceElements_(),
elementMap_(),
comm_(),
pointsRepRef_(),
pointsUniRef_(),
mapUniqueP2Map_(),
mapRepeatedP2Map_(),
AABBTree_(),
volumeID_(10),
edgeElements_(),
surfaceEdgeElements_(),
meshFileName_("fileName.mesh"),
delimiter_(" ")
//#ifdef FULL_TIMER
//,TotalP1ToP2Time_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Total P1 to P2")),
//BuildRedundantInfoTime_ (Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Build Redundant Info")),
//SortUniqueRedundantInfoTime_ (Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Sort Unique redundant")),
//CheckingFlaggedAndDeleteTime_ (Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior")),
//CheckingFlaggedTime_ (Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior: Flagged")),
//DeleteInteriorTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior: Interior")),
//SetGlobalInterfaceIDTime_ (Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Set Global IDs")),
//SetP2InfoTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Set P2 Information")),
//GatherAllMarkedTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior: Flagged: Gather All")),
//UniqueAndFindMarkedTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior: Flagged: Unique and MaxAll")),
//UniqueNewFaceTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior: Flagged: Unique new Faces")),
//SumAllMarkedTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior: Flagged: Sum All")),
//SetP2PointsTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Set P2 Information: Set Points")),
//SetP2ElementsTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Set P2 Information: Set Elements"))
//#endif
{
    edgeElements_ = Teuchos::rcp( new EdgeElements_Type() );
    surfaceEdgeElements_ = Teuchos::rcp( new Elements_Type() );
    surfaceElements_.reset(new Elements());
    elementsC_.reset(new Elements());
    
}

template <class SC, class LO, class GO, class NO>
MeshUnstructured<SC,LO,GO,NO>::MeshUnstructured(CommConstPtr_Type comm, int volumeID):
meshInterface_(),
numElementsGlob_(0),
mapUnique_(),
mapRepeated_(),
pointsRep_(),
pointsUni_(),
bcFlagRep_(),
bcFlagUni_(),
surfaceElements_(),
elementMap_(),
edgeMap_(),
comm_(comm),
pointsRepRef_(),
pointsUniRef_(),
mapUniqueP2Map_(),
mapRepeatedP2Map_(),
AABBTree_(),
volumeID_(volumeID),
edgeElements_(),
surfaceEdgeElements_(),
meshFileName_("fileName.mesh"),
delimiter_(" ")
{
    edgeElements_ = Teuchos::rcp( new EdgeElements_Type() );
    surfaceEdgeElements_ = Teuchos::rcp( new Elements_Type() );
    surfaceElements_.reset(new Elements());
   AABBTree_.reset(new AABBTree_Type());

    elementsC_.reset(new Elements());
        
}

template <class SC, class LO, class GO, class NO>
MeshUnstructured<SC,LO,GO,NO>::~MeshUnstructured(){
}

//template <class SC, class LO, class GO, class NO>
//vec2D_int_ptr_Type MeshUnstructured<SC,LO,GO,NO>::getElements(){
//
//
//    this->elementsVec_ = Teuchos::rcp( new vec2D_int_Type( this->elementsC_->numberElements() ) );
//    for (int i=0; i<this->elementsVec_->size(); i++)
//        this->elementsVec_->at(i) = this->elementsC_->getElement(i).getVectorNodeList();
//
//    return this->elementsVec_;
//}

// ---------------------------------------------------------------------------------------------------
// NEW SETTER AND GETTER
// ---------------------------------------------------------------------------------------------------
template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::setParameterList( ParameterListPtr_Type& pL ) {
    pList_ = pL;
}

template <class SC, class LO, class GO, class NO>
ParameterListConstPtr_Type MeshUnstructured<SC,LO,GO,NO>::getParameterList( ) const{
    return pList_;
}
    
template <class SC, class LO, class GO, class NO>
vec_int_ptr_Type MeshUnstructured<SC,LO,GO,NO>::getElementsFlag() const{
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "we are not using the correct flags here. use the flags of elementC_." );
    vec_int_ptr_Type tmp;
    return tmp;
}

template <class SC, class LO, class GO, class NO>
typename MeshUnstructured<SC,LO,GO,NO>::MapConstPtr_Type MeshUnstructured<SC,LO,GO,NO>::getMapUnique() const{

    return mapUnique_;
}

template <class SC, class LO, class GO, class NO>
typename MeshUnstructured<SC,LO,GO,NO>::MapConstPtr_Type MeshUnstructured<SC,LO,GO,NO>::getMapRepeated() const{

    return mapRepeated_;
}

template <class SC, class LO, class GO, class NO>
typename MeshUnstructured<SC,LO,GO,NO>::MapConstPtr_Type MeshUnstructured<SC,LO,GO,NO>::getMapUniqueP2() const{

    return mapUniqueP2Map_;
}

template <class SC, class LO, class GO, class NO>
typename MeshUnstructured<SC,LO,GO,NO>::MapConstPtr_Type MeshUnstructured<SC,LO,GO,NO>::getMapRepeatedP2() const{

    return mapRepeatedP2Map_;
}

template <class SC, class LO, class GO, class NO>
typename MeshUnstructured<SC,LO,GO,NO>::MapConstPtr_Type MeshUnstructured<SC,LO,GO,NO>::getElementMap(){
    TEUCHOS_TEST_FOR_EXCEPTION( elementMap_.is_null(), std::runtime_error, "Element map of mesh does not exist." );
    return elementMap_;
}

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// edgeMap
template <class SC, class LO, class GO, class NO>
typename MeshUnstructured<SC,LO,GO,NO>::MapConstPtr_Type MeshUnstructured<SC,LO,GO,NO>::getEdgeMap(){
    TEUCHOS_TEST_FOR_EXCEPTION( edgeMap_.is_null(), std::runtime_error, "Edge map of mesh does not exist." );
    return edgeMap_;
}

template <class SC, class LO, class GO, class NO>
vec2D_dbl_ptr_Type MeshUnstructured<SC,LO,GO,NO>::getPointsRepeated() const{

    return pointsRep_;
}

template <class SC, class LO, class GO, class NO>
vec2D_dbl_ptr_Type MeshUnstructured<SC,LO,GO,NO>::getPointsUnique() const{

    return pointsUni_;
}

template <class SC, class LO, class GO, class NO>
vec_int_ptr_Type MeshUnstructured<SC,LO,GO,NO>::getBCFlagRepeated() const{

    return bcFlagRep_;
}

template <class SC, class LO, class GO, class NO>
vec_int_ptr_Type MeshUnstructured<SC,LO,GO,NO>::getBCFlagUnique() const{

    return bcFlagUni_;
}

template <class SC, class LO, class GO, class NO>
typename MeshUnstructured<SC,LO,GO,NO>::ElementsPtr_Type MeshUnstructured<SC,LO,GO,NO>::getElementsC(){
    return elementsC_;
}

template <class SC, class LO, class GO, class NO>
typename MeshUnstructured<SC,LO,GO,NO>::ElementsPtr_Type MeshUnstructured<SC,LO,GO,NO>::getSurfaceElements(){
    return surfaceElements_;
}

template <class SC, class LO, class GO, class NO>
int MeshUnstructured<SC,LO,GO,NO>::getDimension(){

    return dim_;
}

template <class SC, class LO, class GO, class NO>
GO MeshUnstructured<SC,LO,GO,NO>::getNumElementsGlobal(){

    return numElementsGlob_;
}

template <class SC, class LO, class GO, class NO>
LO MeshUnstructured<SC,LO,GO,NO>::getNumElements(){
    TEUCHOS_TEST_FOR_EXCEPTION( this->elementsC_.is_null(), std::runtime_error ,"Elements do not exist." );
    return this->elementsC_->numberElements();
}

template <class SC, class LO, class GO, class NO>
LO MeshUnstructured<SC,LO,GO,NO>::getNumPoints(std::string type){
    if (!type.compare("Unique"))
    return pointsUni_->size();
    else if(!type.compare("Repeated"))
    return pointsRep_->size();

    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Select valid map type: unique or repeated.");
    return 0;
}

template <class SC, class LO, class GO, class NO>
int MeshUnstructured<SC,LO,GO,NO>::getOrderElement(){

    switch (dim_) {
        case 2:
            if ( !FEType_.compare("P1") )
                return 3;
            else if ( !FEType_.compare("P1-disc") || !FEType_.compare("P1-disc-global") ){
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "P1-disc only available in 3D.");
            }
            else if( !FEType_.compare("P2") )
                return 6;
            else if( !FEType_.compare("P2-CR") ){
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "P2-CR only available in 3D.");
            }
            else if( !FEType_.compare("Q2-20") ){
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Q2-20 only available in 3D.");
            }
            else if( !FEType_.compare("Q2") ){
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Q2 only available in 3D.");
            }
            break;
        case 3:
            if ( !FEType_.compare("P1") )
                return 4;
            if ( !FEType_.compare("P1-disc") || !FEType_.compare("P1-disc-global") )
                return 4;
            else if( !FEType_.compare("P2") )
                return 10;
            else if( !FEType_.compare("P2-CR") )
                return 15;
            else if( !FEType_.compare("Q2-20") )
                return 20;
            else if( !FEType_.compare("Q2") )
                return 27;
            break;
        default:
            return -1;
            break;
    }
    return -1;
}



template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::getTriangles(int vertex1ID, int vertex2ID, vec_int_Type &vertices3ID){
    vertices3ID.resize(2);
    if (vertex2ID<vertex1ID) {
        int tmp = vertex2ID;
        vertex2ID = vertex1ID;
        vertex1ID = tmp;
    }
    if (vertex1ID == 0) {
        if (vertex2ID == 1) {
            vertices3ID.at(0) = 2;
            vertices3ID.at(1) = 3;
        }
        else if (vertex2ID == 2) {
            vertices3ID.at(0) = 1;
            vertices3ID.at(1) = 3;
        }
        else if (vertex2ID == 3) {
            vertices3ID.at(0) = 1;
            vertices3ID.at(1) = 2;
        }

    }
    else if(vertex1ID == 1){
        if (vertex2ID == 2) {
            vertices3ID.at(0) = 0;
            vertices3ID.at(1) = 3;
        }
        else if (vertex2ID == 3) {
            vertices3ID.at(0) = 0;
            vertices3ID.at(1) = 2;
        }
    }

    else if(vertex1ID == 2){
        if (vertex2ID == 3) {
            vertices3ID.at(0) = 0;
            vertices3ID.at(1) = 1;
        }
    }
    else {
#ifdef ASSERTS_WARNINGS
        MYASSERT(false,"Could not identify triangle.");
#endif
    }
}


/*template <class SC, class LO, class GO, class NO>
int MeshUnstructured<SC,LO,GO,NO>::determineFlagP2( FiniteElement& fe, LO p1ID, LO p2ID, vec2D_int_Type& permutation ){

    int newFlag = std::numeric_limits<int>::max();
    vec_int_Type newFlags(0);
    bool foundLineSegment;
    const vec_int_Type nodeList = fe.getVectorNodeList();
    
    vec_int_Type localElementNumbering(2);
    auto it1 = find( nodeList.begin(), nodeList.end() , p1ID );
    localElementNumbering[0] = distance( nodeList.begin() , it1 );
    auto it2 = find( nodeList.begin(), nodeList.end() , p2ID );
    localElementNumbering[1] = distance( nodeList.begin() , it2 );

    fe.findEdgeFlagInSubElements( localElementNumbering, newFlags, false , permutation, foundLineSegment );
    
    if (newFlags.size() == 0)
        newFlag = this->volumeID_;
    
    // We use the lowest flag
    for (int k = 0; k < newFlags.size(); k++) {
        if ( newFlag > newFlags[k] )
            newFlag = newFlags[k];
    }
    
    return newFlag;
}*/
    
template <class SC, class LO, class GO, class NO>
int MeshUnstructured<SC,LO,GO,NO>::determineFlagP2( LO p1ID, LO p2ID, LO localEdgeID, vec2D_LO_Type& markedPoint ){
    

    ElementsPtr_Type elements = this->getElementsC();

    vec2D_int_Type permutation = elements->getElementEdgePermutation();

    EdgeElementsPtr_Type edgeElements = this->getEdgeElements();

    MapConstPtr_Type elementMap = this->getElementMap();

    int rank = this->comm_->getRank();

    int flag1 = ( *this->getBCFlagRepeated() )[p1ID];

    int flag2 = ( *this->getBCFlagRepeated() )[p2ID];

    int newFlag = std::numeric_limits<int>::max();

    if(flag1 == this->volumeID_ || flag2 == this->volumeID_ ) // one node is in an inner node, than the new node is an inner node aswell
        newFlag = this->volumeID_;
    else{

        // check if node 1 and node 2 are part of the same surface. In this case we can use the flag of the corresponding surface. Otherwise we mark the new point as an interior point and it gets the volumeID_.

        const vec_LO_Type elementsOfEdge = edgeElements->getElementsOfEdge( (int) localEdgeID );

        const vec_GO_Type elementsOfEdgeGlobal = edgeElements->getElementsOfEdgeGlobal( (int) localEdgeID );


        bool markPoint = false;
        bool foundFlag = false;
        vec_int_Type localElementNumbering(2);
        vec_int_Type edge(2);
        bool foundLineSegment = false;
        vec_int_Type newFlags(0);
        for (int i=0; i<elementsOfEdge.size() && !foundLineSegment; i++) {
            if ( elementsOfEdge[i] != OTLO::invalid() ) {
                //In elementsOfEdge we can access all elements which have this (the current) edge.
                
                FiniteElement fe = elements->getElement( elementsOfEdge[i] );
                const vec_int_Type nodeList = fe.getVectorNodeList();

                // we need to determine the numbering of p1ID and p2ID corresponding to the current element
                auto it1 = find( nodeList.begin(), nodeList.end() , p1ID );
                localElementNumbering[0] = distance( nodeList.begin() , it1 );
                auto it2 = find( nodeList.begin(), nodeList.end() , p2ID );
                localElementNumbering[1] = distance( nodeList.begin() , it2 );
                edge[0] = p1ID; edge[1] = p2ID;
                sort( edge.begin(), edge.end() );
                
                fe.findEdgeFlagInSubElements( edge, newFlags, false /*we are not in a subElement yet*/, permutation, foundLineSegment );

                //We need to mark this point since it can still be on the surface and another element holds the corresponding surface with the correct flag.
                if (newFlags.size() == 0 && newFlag > this->volumeID_)
                    newFlag = this->volumeID_; //do we need this?

                //If we found a line element, the we choose this flag
                if (foundLineSegment){
                    foundFlag = true;
                    newFlag = newFlags [0];
                }
                else {
                    // We use the lowest flag of all surfaces
                    for (int k = 0; k < newFlags.size(); k++) {
                        foundFlag = true;
                        if (newFlag > newFlags[k] )
                            newFlag = newFlags[k];
                    }
                }
            }
            else
                markPoint = true;
        }
        if (markPoint && !foundFlag) {
            // We have not found a surface flag for this point.
            // No part of the element which owns this edge is part of the surface, but there might be a different element which has this edge but does not own it, but this element might have a part of the surface. We mark this point and determine the correct flag later
            markedPoint.push_back( {p1ID, p2ID, localEdgeID } );
            newFlag = -1;
        }
        
    }
    
    return newFlag;
}

// ------------------------


template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::findSurfaces( const vec_int_Type& elementNodeList, vec_int_Type numbering,  vec2D_int_Type& localSurfaceNodeList_vec, vec_int_Type& loc_vector, bool critical){
    int tmpID;
    loc_vector.resize(0);
    vec2D_int_Type::iterator it;
    int id1 = elementNodeList.at(numbering.at(0));
    int id2 = elementNodeList.at(numbering.at(1));

    if (this->dim_==2) {
        vec_int_Type searchfor(2,-1);
        // Here, we use that the local surfaces are sorted. See MyMeshConvert::ReadMesh(...)
        if ( id1 < id2 ) {
            searchfor.at(0) = id1;
            searchfor.at(1) = id2;
        }
        else{
            searchfor.at(0) = id2;
            searchfor.at(1) = id1;
        }

        it = find( localSurfaceNodeList_vec.begin(), localSurfaceNodeList_vec.end(), searchfor );
        if (it!=localSurfaceNodeList_vec.end()) {
            loc_vector.push_back( distance(localSurfaceNodeList_vec.begin(),it) );
        }

    }
    else if (this->dim_==3){
        
        vec_int_Type vertices3ID(0);
        //we want to identify the possible third surface component of a triangle. If we use other elements in 3D we might change this function.
        getTriangles(numbering.at(0), numbering.at(1), vertices3ID);
        // we have the two  possible combinations as a vector: vertices3ID
        
        vec_int_Type searchfor(3,-1);
        for (int i=0; i<2; i++) {
            searchfor.at(0) = id1;
            searchfor.at(1) = id2;
            searchfor.at(2) = elementNodeList.at( vertices3ID.at(i) );

            sort( searchfor.begin(), searchfor.end() );
            
            
            it = find( localSurfaceNodeList_vec.begin(), localSurfaceNodeList_vec.end(), searchfor );
            if (it!=localSurfaceNodeList_vec.end()){
                loc_vector.push_back( distance(localSurfaceNodeList_vec.begin(),it) );
            }
        }
    }
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::findEdges( const vec_int_Type& elementNodeList, vec_int_Type numbering,  vec2D_int_Type& localEdgeNodeList_vec, vec_int_Type& loc_vector){
    int tmpID;
    loc_vector.resize(0);
    vec2D_int_Type::iterator it;
    int id1 = elementNodeList.at(numbering.at(0));
    int id2 = elementNodeList.at(numbering.at(1));
    vec_int_Type searchEdge(0);
    if (id2 > id1)
        searchEdge = { id1, id2 };
    else
        searchEdge = { id2, id1 };
    
    it = find( localEdgeNodeList_vec.begin(), localEdgeNodeList_vec.end(), searchEdge );
    if ( it != localEdgeNodeList_vec.end() )
        loc_vector.push_back( distance(localEdgeNodeList_vec.begin(),it) );
}

template <class SC, class LO, class GO, class NO>
typename MeshUnstructured<SC,LO,GO,NO>::MeshInterfacePtr_Type MeshUnstructured<SC,LO,GO,NO>::getMeshInterface(){
    TEUCHOS_TEST_FOR_EXCEPTION( meshInterface_.is_null(), std::runtime_error, "MeshInterface is null.");
    return meshInterface_;
}
    
template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::setMeshFileName(string meshFileName, string delimiter){

    meshFileName_ = meshFileName;
    delimiter_ = delimiter;
}

// Mesh Functions
template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::setElementFlags(std::string type){

    ElementsPtr_Type elements = this->getElementsC();
//    this->elementFlag_.reset( new vec_int_Type( elements->numberElements(), 0 ) );
    if (type == "TPM_square") {
        double xRef, yRef;

        for (int i=0; i<elements->numberElements(); i++) {
            xRef = ( this->pointsRep_->at( elements->getElement(i).getNode(0) )[0] + this->pointsRep_->at( elements->getElement(i).getNode(1) )[0] + this->pointsRep_->at( elements->getElement(i).getNode(2) )[0] ) / 3.;
            yRef = ( this->pointsRep_->at( elements->getElement(i).getNode(0) )[1] + this->pointsRep_->at( elements->getElement(i).getNode(1) )[1] + this->pointsRep_->at( elements->getElement(i).getNode(2) )[1] ) / 3.;
            if ( xRef>=0.3  && xRef<=0.7) {
                if ( yRef>= 0.6) {
                    elements->getElement(i).setFlag(1);
//                    this->elementFlag_->at(i) = 1;
                }
            }
        }
    }
    else if (type == "Excavation1"){
    }
    else{
    
    }

}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::setReferenceConfiguration()
{
    // Bemerkung: Repeated und Unique sind unterschiedlich lang!!! => zwei Schleifen

    // Setze zunaechst alles auf Null, andernfalls kann man nicht drauf zugreifen
    //    vec2D_dbl_ptr_Type zeroRep(new vec2D_dbl_Type(pointsRep_->size(),vec_dbl_Type(pointsRep_->at(0).size(),0.0)));
    //    vec2D_dbl_ptr_Type zeroUni(new vec2D_dbl_Type(pointsUni_->size(),vec_dbl_Type(pointsUni_->at(0).size(),0.0)));


    pointsRepRef_.reset( new vec2D_dbl_Type() );
    pointsUniRef_.reset( new vec2D_dbl_Type() );
    // zeroRep und zeroUni leben nur hier drinnen, weswegen wir Pointer gleich Pointer setzen koennen.
    // TODO: *PointsRepRef_ = *zeroRep funktioniert nicht.
    *pointsRepRef_ = *pointsRep_;
    *pointsUniRef_ = *pointsUni_;

    //    // Repeated
    //    for(int i = 0; i < pointsRep_->size(); i++)
    //    {
    //        for(int j = 0; j < pointsRep_->at(0).size(); j++)
    //        {
    //            pointsRepRef_->at(i).at(j) = PointsRep_->at(i).at(j);
    //        }
    //    }
    //
    //    // Unique
    //    for(int i = 0; i < PointsUni_->size(); i++)
    //    {
    //        for(int j = 0; j < PointsUni_->at(0).size(); j++)
    //        {
    //            PointsUniRef_->at(i).at(j) = PointsUni_->at(i).at(j);
    //        }
    //    }
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::create_AABBTree(){
    if (AABBTree_.is_null()){
        AABBTree_.reset(new AABBTree_Type());
    }
    AABBTree_->createTreeFromElements(
        getElementsC(),
        getPointsRepeated()
    );
}

template <class SC, class LO, class GO, class NO>
vec_int_ptr_Type MeshUnstructured<SC,LO,GO,NO>::findElemsForPoints(
    vec2D_dbl_ptr_Type queryPoints
){
    int numPoints = queryPoints->size();

    // Return vector. -1 means that point is in no elem, otherwise entry is the
    // elem the point is in
    vec_int_ptr_Type pointToElem(
        new vec_int_Type(
            numPoints,
            -1
        )
    );

    // Create tree if it is empty
    if (AABBTree_->isEmpty()){
        AABBTree_->createTreeFromElements(
            getElementsC(),
            getPointsRepeated(),
            false
        );
    }

    // Query the AABBTree
    map<int, list<int> > treeToItem;
    map<int, list<int> > itemToTree;
    tie(treeToItem, itemToTree) = AABBTree_->scanTree(queryPoints, false);

    // FIXME: put this in a function of AABBTree?
    // unnest the returned answer for each query_point
    int point = -1;
    bool found = false;
    list<int> rectangles;
    list<int> elements;
    for (auto keyValue: itemToTree){
        // FIXME: put this in a function of AABBTree?
        // rectangles is a list<int> of all rectangles point is in
        // find the element(s) that is/are in all of said rectangles.
        // If there is only one element that is the element the point is in,
        // if not we have to query all remaining elements
        point = keyValue.first;
        rectangles = keyValue.second;


        // query all remaining elements
        for(auto rectangle: rectangles){
            elements = AABBTree_->getElements(rectangle);
            for(auto element: elements){
                found = isPointInElem(queryPoints->at(point), element);
                if( found ){
                    pointToElem->at(point) = element;
                    break;
                }
            }
            if( found ){
                // we already found the element, no need to check additional rectangles
                break;
            }
        }

    }
    return pointToElem;
}

template <class SC, class LO, class GO, class NO>
vec_dbl_Type MeshUnstructured<SC,LO,GO,NO>::getBaryCoords(vec_dbl_Type point, int element){
    vec_int_Type localNodes = elementsC_->getElement(element).getVectorNodeList();
    vec2D_dbl_Type coords(
        localNodes.size(),
        vec_dbl_Type(2, 0.0) // FIXME: this depends on the dimension
    );
    int node = 0;
    for (int localNode=0; localNode<localNodes.size(); localNode++){
        node = localNodes.at(localNode); // get global node_id
        coords.at(localNode) = pointsRep_->at(node); //get global coordinates
    }

    double px, py, x1, x2, x3, y1, y2, y3;
    px = point.at(0); py = point.at(1);
    x1 = coords.at(0).at(0); y1 = coords.at(0).at(1);
    x2 = coords.at(1).at(0); y2 = coords.at(1).at(1);
    x3 = coords.at(2).at(0); y3 = coords.at(2).at(1);

    // baryzentric coordinates
    double det_T = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);

    vec_dbl_Type baryCoords(
        3,
        0.0
    );
    baryCoords[0] =(y2 - y3) * (px - x3) + (x3 - x2) * (py - y3);
    baryCoords[0] = baryCoords[0] / det_T;

    baryCoords[1] = (y3 -y1) * (px - x3) + (x1 - x3) * (py - y3);
    baryCoords[1] = baryCoords[1] / det_T;

    baryCoords[2] = 1 - baryCoords[1] - baryCoords[0];
    return baryCoords;
}


template <class SC, class LO, class GO, class NO>
bool MeshUnstructured<SC,LO,GO,NO>:: isPointInElem(vec_dbl_Type point, int element){
    // FIXME: This is only valid for a triangle
    vec_dbl_Type baryCoords;
    baryCoords = getBaryCoords(point, element);

    if(baryCoords[0] >= 0 && baryCoords[1] >= 0 && baryCoords[2] >= 0){
        return true;
    }
    return false;
}
/*!

 \brief Building edgeMap after refinement.

@param[in] mapGlobalProc Map of global processor numbers
@param[in] mapProc Map of local processor number

*/
template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::buildEdgeMap(){

		int maxRank = std::get<1>(this->rankRange_);
		vec_GO_Type globalProcs(0);
		for (int i=0; i<= maxRank; i++)
			globalProcs.push_back(i);

		Teuchos::ArrayView<GO> globalProcArray = Teuchos::arrayViewFromVector( globalProcs);

		vec_GO_Type localProc(0);
		localProc.push_back(this->comm_->getRank());
		Teuchos::ArrayView<GO> localProcArray = Teuchos::arrayViewFromVector( localProc);

		MapPtr_Type mapGlobalProc =
			Teuchos::rcp( new Map_Type( this->getEdgeMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalProcArray, 0, this->comm_) );

		MapPtr_Type mapProc =
			Teuchos::rcp( new Map_Type( this->getEdgeMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), localProcArray, 0, this->comm_) );
		

		vec2D_int_Type interfaceEdgesLocalId(1,vec_int_Type(1));
		const int myRank = this->comm_->getRank();

		MultiVectorLOPtr_Type exportLocalEntry = Teuchos::rcp( new MultiVectorLO_Type( mapProc, 1 ) );

		// (A) First we determine a Map only for the interface Nodes
		// This will reduce the size of the Matrix we build later significantly if only look at the interface edges
		int numEdges= this->edgeElements_->numberElements();
		vec2D_GO_Type inzidenzIndices(0,vec_GO_Type(2)); // Vector that stores global IDs of each edge (in Repeated Sense)
		vec_LO_Type localEdgeIndex(0); // stores the local ID of edges in question 
		vec_GO_Type id(2);
		int edgesUnique=0;
    	EdgeElementsPtr_Type edgeElements = this->edgeElements_; // Edges

		vec2D_dbl_ptr_Type points = this->pointsRep_;

		int interfaceNum=0;
		for(int i=0; i<numEdges; i++ ){
			if(edgeElements->getElement(i).isInterfaceElement()){

				id[0] = this->mapRepeated_->getGlobalElement(edgeElements->getElement(i).getNode(0)); 
				id[1] = this->mapRepeated_->getGlobalElement(edgeElements->getElement(i).getNode(1));
			 	


				sort(id.begin(),id.end());
				inzidenzIndices.push_back(id);

				localEdgeIndex.push_back(i);
				interfaceNum++;
			}
	
			else{
				edgesUnique++;
			}


		 }
		// This Matrix is row based, where the row is based on mapInterfaceNodesUnqiue
		// We then add a '1' Entry when two global Node IDs form an edge
		MatrixPtr_Type inzidenzMatrix = Teuchos::rcp( new Matrix_Type(this->mapUnique_, 40 ) );
		Teuchos::Array<GO> index(1);
		Teuchos::Array<GO> col(1);
		Teuchos::Array<SC> value(1, Teuchos::ScalarTraits<SC>::one() );

		for(int i=0; i<inzidenzIndices.size(); i++ ){
	
			index[0] = inzidenzIndices[i][0];
			col[0] = inzidenzIndices[i][1];
			inzidenzMatrix->insertGlobalValues(index[0], col(), value());
		
		 }
   		inzidenzMatrix->fillComplete(); //mapInterfaceNodesUnique,mapInterfaceNodesUnique);
		
	
		// ---------------------------------------------------
		// 2 .Set unique edges IDs ---------------------------
		// Setting the IDs of Edges that are uniquely on one
		// Processor
		// ---------------------------------------------------
		exportLocalEntry->putScalar( (LO) edgesUnique );

		MultiVectorLOPtr_Type newEdgesUniqueGlobal= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		newEdgesUniqueGlobal->putScalar( (LO) 0 ); 
		newEdgesUniqueGlobal->importFromVector( exportLocalEntry, true, "Insert");
		// offset EdgesUnique for proc and globally
		Teuchos::ArrayRCP< const LO > newEdgesList = newEdgesUniqueGlobal->getData(0);

		GO procOffsetEdges=0;
		for(int i=0; i< myRank; i++)
			procOffsetEdges= procOffsetEdges + newEdgesList[i];

		// global IDs for map
		vec_GO_Type vecGlobalIDsEdges(this->edgeElements_->numberElements()); 
	
		// Step 1: adding unique global edge IDs
		int count=0;
		for(int i=0; i< this->edgeElements_->numberElements(); i++){
			if(!this->edgeElements_->getElement(i).isInterfaceElement()){
				vecGlobalIDsEdges.at(i) = procOffsetEdges+count;
				count++;
			}
		}	
		
		// Now we add the repeated ids, by first turning interfaceEdgesTag into a map
		// Offset for interface IDS:
		GO offsetInterface=0;
		for(int i=0; i< maxRank+1; i++)
			 offsetInterface=  offsetInterface + newEdgesList[i];
		
		//Now we count the row entries on each processor an set global IDs

		Teuchos::ArrayView<const LO> indices;
		Teuchos::ArrayView<const SC> values;
		vec2D_GO_Type inzidenzIndicesUnique(0,vec_GO_Type(2)); // Vector that stores only both global IDs if the first is part of my unique Interface Nodes
		MapConstPtr_Type colMap = inzidenzMatrix->getMap("col");
		MapConstPtr_Type rowMap = inzidenzMatrix->getMap("row");
		int numRows = rowMap->getNodeNumElements();
		int uniqueEdges =0;
		for(int i=0; i<numRows; i++ ){
			inzidenzMatrix->getLocalRowView(i, indices,values); 
			uniqueEdges = uniqueEdges+indices.size();
			vec_GO_Type edgeTmp(2);
			for(int j=0; j<indices.size(); j++){
				edgeTmp[0] = rowMap->getGlobalElement(i);
				edgeTmp[1] = colMap->getGlobalElement(indices[j]);
				inzidenzIndicesUnique.push_back(edgeTmp);
			}
		}
	
		exportLocalEntry->putScalar( uniqueEdges );
		MultiVectorLOPtr_Type newEdgesInterfaceGlobal= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		newEdgesInterfaceGlobal->putScalar( (LO) 0 ); 
		newEdgesInterfaceGlobal->importFromVector( exportLocalEntry, true, "Insert");

		// offset EdgesUnique for proc and globally
		Teuchos::ArrayRCP< const LO > numUniqueInterface = newEdgesInterfaceGlobal->getData(0);

		procOffsetEdges=0;
		for(int i=0; i< myRank; i++)
			procOffsetEdges= procOffsetEdges + numUniqueInterface[i];

		int numInterfaceEdges=0;
		
		vec_GO_Type uniqueInterfaceIDsList_(inzidenzIndicesUnique.size());
		for(int i=0; i< uniqueInterfaceIDsList_.size(); i++)
			uniqueInterfaceIDsList_[i] = procOffsetEdges + i;

		MatrixPtr_Type indMatrix = Teuchos::rcp( new Matrix_Type(this->mapUnique_, 40 ) );

		for(int i=0; i<inzidenzIndicesUnique.size(); i++ ){
			index[0] = inzidenzIndicesUnique[i][0];
			col[0] = inzidenzIndicesUnique[i][1];
			Teuchos::Array<SC> value2(1,uniqueInterfaceIDsList_[i]);
			indMatrix->insertGlobalValues(index[0], col(), value2());
		 }
   		indMatrix->fillComplete(); 

		MatrixPtr_Type importMatrix = Teuchos::rcp( new Matrix_Type(this->mapRepeated_, 40 ) );
   		
		importMatrix->importFromVector(indMatrix,false,"Insert");
		importMatrix->fillComplete(); 		
		
		// Determine global indices
		GO edgeID=0;
		colMap = importMatrix->getMap("col");
		rowMap = importMatrix->getMap("row");
	
		LO valueID=0;
		bool found = false;
		GO entry =0;
		for(int i=0; i<inzidenzIndices.size(); i++ ){
			
			importMatrix->getLocalRowView(rowMap->getLocalElement(inzidenzIndices[i][0]), indices,values); // Indices and values connected to node i / row i in Matrix
			// Entries in 'indices' represent the local entry in 'colmap
			// with 'getGlobalElement' we know the global Node ID that belongs to the first Node that form an edge
			// vector in with entries only for edges belonging to node i;
			vec2D_GO_Type indicesTmp(indices.size(),vec_GO_Type(2));
			vec_GO_Type indTmp(2);
			for(int j=0; j<indices.size(); j++){
				indTmp[0] = colMap->getGlobalElement(indices[j]);
				indTmp[1] = values[j];
				indicesTmp.push_back(indTmp);	// vector with the indices and values belonging to node i
			}
			//sort(indicesTmp.begin(),indicesTmp.end());
			found = false;
			for(int k=0; k<indicesTmp.size();k++){
				if(inzidenzIndices[i][1] == indicesTmp[k][0]){
					entry =k;
					k = indicesTmp.size();
					edgeID = indicesTmp[entry][1];
					vecGlobalIDsEdges.at(localEdgeIndex[i]) = offsetInterface + edgeID;
					found =true;
				}
			}
			if(found == false)
				cout << " Asking for row " << rowMap->getLocalElement(inzidenzIndices[i][0]) << " for Edge [" << inzidenzIndices[i][0] << ",  " << inzidenzIndices[i][1] << "], on Proc " << myRank << " but no Value found " <<endl;
		 }


		Teuchos::RCP<std::vector<GO>> edgesGlobMapping = Teuchos::rcp( new vector<GO>( vecGlobalIDsEdges ) );
		Teuchos::ArrayView<GO> edgesGlobMappingArray = Teuchos::arrayViewFromVector( *edgesGlobMapping);

		this->edgeMap_.reset(new Map<LO,GO,NO>(this->getMapRepeated()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesGlobMappingArray, 0, this->comm_) );
		//this->edgeMap_->print();
}
/*!

 \brief Updating ElementsOfEdgesLocal and ElementsOfEdgesGlobal.

@param[in] maxRank The maximal processor rank.
@param[in] edgeMap Map of global edge ids.

*/

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::updateElementsOfEdgesLocalAndGlobal(){

	int maxRank = std::get<1>(this->rankRange_);
	if(maxRank >0 && this->dim_ == 2){
		vec_GO_Type edgesInterfaceGlobalID(0);
		LO id=0;
		for(int i=0; i< this->edgeElements_->numberElements(); i++){
			if(this->edgeElements_->getElement(i).isInterfaceElement() ){		
				this->edgeElements_->setElementsOfEdgeLocalEntry(i,-1);
				edgesInterfaceGlobalID.push_back(this->edgeMap_->getGlobalElement(i)); // extracting the global IDs of the new interfaceEdges
			}
			
		}

		// communticating elements across interface
		Teuchos::ArrayView<GO> edgesInterfaceGlobalID_ = Teuchos::arrayViewFromVector( edgesInterfaceGlobalID);

		MapPtr_Type mapGlobalInterface =
			Teuchos::rcp( new Map_Type( this->edgeMap_->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesInterfaceGlobalID_, 0, this->comm_) );
		//mapGlobalInterface->print();

		// Global IDs of Procs
		// Setting newPoints as to be communicated Values
		MultiVectorLOPtr_Type interfaceElements = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		Teuchos::ArrayRCP< LO > interfaceElementsEntries  = interfaceElements->getDataNonConst(0);

		for(int i=0; i< interfaceElementsEntries.size() ; i++){
			interfaceElementsEntries[i] = this->edgeElements_->getElementsOfEdgeGlobal(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i])).at(0);
		}

		MapConstPtr_Type mapGlobalInterfaceUnique = mapGlobalInterface->buildUniqueMap( this->rankRange_ );

		MultiVectorLOPtr_Type isInterfaceElement_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_imp->putScalar( (LO) 0 ); 
		isInterfaceElement_imp->importFromVector( interfaceElements, false, "Insert");

		MultiVectorLOPtr_Type isInterfaceElement_exp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_exp->putScalar( (LO) 0 ); 
		isInterfaceElement_exp->exportFromVector( interfaceElements, false, "Insert");

		MultiVectorLOPtr_Type isInterfaceElement2_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		isInterfaceElement2_imp->putScalar( (LO) 0 ); 
		isInterfaceElement2_imp->importFromVector(isInterfaceElement_imp, false, "Insert");

		isInterfaceElement2_imp->exportFromVector(isInterfaceElement_exp, false, "Insert");

		interfaceElementsEntries  = isInterfaceElement2_imp->getDataNonConst(0);

		for(int i=0; i< interfaceElementsEntries.size() ; i++){
			this->edgeElements_->setElementsOfEdgeGlobalEntry(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i]),interfaceElementsEntries[i]);
			}

	}

	// Contrary to the 2D case in 3D an edge can be part of more than two elements
	// We need to determine how many elements are connected to an edge and import the global IDs from different Processors
	if(maxRank >0 && this->dim_ == 3){

		// First we determine the interface edges, which entries in elementsOfEdgesGlobal/Local have to be completed
		vec_GO_Type edgesInterfaceGlobalID(0);
		vec_int_Type numberElements(0); // represents number of elements the edge is connected to on my Processor
		
		for(int i=0; i< this->edgeElements_->numberElements(); i++){
			if(this->edgeElements_->getElement(i).isInterfaceElement() ){
				edgesInterfaceGlobalID.push_back(this->edgeMap_->getGlobalElement(i)); // extracting the global IDs of the new interfaceEdges
			}	
		}
		sort(edgesInterfaceGlobalID.begin(), edgesInterfaceGlobalID.end());

		for(int i=0; i< edgesInterfaceGlobalID.size(); i++){
			numberElements.push_back(this->edgeElements_->getElementsOfEdgeGlobal(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i])).size());
		}
			
		
		// from this we build a map
		Teuchos::ArrayView<GO> edgesInterfaceGlobalID_ = Teuchos::arrayViewFromVector( edgesInterfaceGlobalID);

		MapPtr_Type mapGlobalInterface =
			Teuchos::rcp( new Map_Type( this->edgeMap_->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesInterfaceGlobalID_, 0, this->comm_) );

		// As edges can be part of multiple elements on different processors we collect the number of elements connected to the edge in total
		MultiVectorLOPtr_Type numberInterfaceElements = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		Teuchos::ArrayRCP< LO > numberInterfaceElementsEntries  = numberInterfaceElements->getDataNonConst(0);

		for(int i=0; i< numberInterfaceElementsEntries.size(); i++)
			numberInterfaceElementsEntries[i] = numberElements[i];

		MapConstPtr_Type mapGlobalInterfaceUnique = mapGlobalInterface->buildUniqueMap( this->rankRange_ );
	
		// Element are unique to each processor. This means that the number we collect is the number of elements that are connected to my edge on other processors.
		// With the following communication we add up all the entries for a certain global Edge ID
		// Potential causes of error:
		//  - if an edge is not identified as an interface edge, of course it will not import nor export its interface Information, making itself and others incomplete

		MultiVectorLOPtr_Type isInterfaceElement_exp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_exp->putScalar( (LO) 0 ); 
		isInterfaceElement_exp->exportFromVector( numberInterfaceElements, false, "Add");

		MultiVectorLOPtr_Type isInterfaceElement2_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		isInterfaceElement2_imp->putScalar( (LO) 0 ); 
		isInterfaceElement2_imp->importFromVector(isInterfaceElement_exp, true, "Insert");

		Teuchos::ArrayRCP< LO > numberInterfaceElementsImportEntries  = isInterfaceElement2_imp->getDataNonConst(0);

		vec_int_Type missingEntries(numberInterfaceElementsEntries.size());
		// With this number we can complete the elementsOfEdgeLocal List with -1 for the elements not on our processor
		for(int i=0; i<numberInterfaceElementsEntries.size() ; i++){
			for(int j=0; j< numberInterfaceElementsImportEntries[i] - numberInterfaceElementsEntries[i];j++){
				this->edgeElements_->setElementsOfEdgeLocalEntry(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i]),-1);
				missingEntries[i] = numberInterfaceElementsImportEntries[i] -numberInterfaceElementsEntries[i];
			}
		}

		// Next we need to identify the global Element IDs of those missing entries and communicate them
		// Hey i got the global Elements ... of edge ... -> exchange
		// Elements are uniquely distributed -> you cannot import an element you already have
		// I need x number of entries -> import all i need, export all i have 
		// Global IDs of Procs
		// Setting newPoints as to be communicated Values

		// Communicating max number of necessary values:
		vec_int_Type::iterator it;
		it = max_element(numberElements.begin(), numberElements.end());
		int myNumberElementsMax = numberElements.at(distance(numberElements.begin(), it)); // accumulate(errorElement.begin(), errorElement.end(),0.0);

		reduceAll<int, int> (*this->comm_, REDUCE_MAX,  myNumberElementsMax , outArg ( myNumberElementsMax));

		MultiVectorLOPtr_Type interfaceElements = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		Teuchos::ArrayRCP< LO > interfaceElementsEntries  = interfaceElements->getDataNonConst(0);

		vec2D_int_Type importElements(this->edgeElements_->getElementsOfEdgeGlobal().size(),vec_int_Type( 0));

		// We extended this function to also consider the ranks. Before we would only exchange the information, but in case one processor received information from more than one processor at 
		// the same time some of the information would get lost. Now we only send the information one processor holds to the other at the same time and move through the processor only destributing their
		// edgeOfElementGlobal Information in a circle like order
		for(int k=0; k< maxRank+1 ; k++){
			
			vec_GO_Type edgesInterfaceGlobalIDProc;				
			if(this->comm_->getRank() == k ){
				edgesInterfaceGlobalIDProc = edgesInterfaceGlobalID; // extracting the global IDs of the new interfaceEdges
			}
					

			// from this we build a map
			Teuchos::ArrayView<GO> edgesInterfaceGlobalIDProc_ = Teuchos::arrayViewFromVector( edgesInterfaceGlobalIDProc);

			MapPtr_Type mapGlobalInterfaceProcs =
				Teuchos::rcp( new Map_Type( this->edgeMap_->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesInterfaceGlobalIDProc_, 0, this->comm_) );

			for(int j=0; j< myNumberElementsMax; j++){
				MultiVectorLOPtr_Type interfaceElements = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceProcs, 1 ) );
				Teuchos::ArrayRCP< LO > interfaceElementsEntries  = interfaceElements->getDataNonConst(0);

				for(int i=0; i< interfaceElementsEntries.size() ; i++){		
					if(numberElements[i] > j && this->comm_->getRank() == k )
						interfaceElementsEntries[i] = this->edgeElements_->getElementsOfEdgeGlobal(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i])).at(j);
					else
						interfaceElementsEntries[i] = -1; 
				}

				MultiVectorLOPtr_Type isInterfaceElement_exp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
				isInterfaceElement_exp->putScalar( (LO) -1 ); 
				isInterfaceElement_exp->exportFromVector( interfaceElements, false, "Insert");

				if(this->comm_->getRank() == k && mapGlobalInterfaceUnique->getNodeNumElements() > 0){
					Teuchos::ArrayRCP< LO > interfaceElementsEntries_exp  = isInterfaceElement_exp->getDataNonConst(0);
					for(int i=0; i<  interfaceElementsEntries_exp.size() ; i++){
						LO id = mapGlobalInterface->getLocalElement(mapGlobalInterfaceUnique->getGlobalElement(i));
						interfaceElementsEntries_exp[i] = interfaceElementsEntries[id];
					}
									
				}

			
				MultiVectorLOPtr_Type isInterfaceElement2_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
				isInterfaceElement2_imp->putScalar( (LO) 0 ); 
				isInterfaceElement2_imp->importFromVector(isInterfaceElement_exp, false, "Insert");

				interfaceElementsEntries  = isInterfaceElement2_imp->getDataNonConst(0);

				for(int i=0; i< interfaceElementsEntries.size() ; i++){
					if(this->comm_->getRank() != k && interfaceElementsEntries[i] != -1)
						importElements[i].push_back( interfaceElementsEntries[i]);
				}
				
			}

		}
		for(int i=0; i< interfaceElementsEntries.size() ; i++){
			sort(importElements[i].begin(),importElements[i].end());
			importElements[i].erase( unique(importElements[i].begin(), importElements[i].end() ), importElements[i].end() );
			if(importElements[i].size() != missingEntries[i])
			   cout << " On Processor " << this->comm_->getRank() << " uneven number for edge imported: " << importElements[i].size() << " missing " << missingEntries[i] << " " << edgesInterfaceGlobalID[i] << endl; // " something went wrong while updating elementsOfEdgesGlobal as the imported entries do not match the supposed number of imported entries. Please check." << endl; 	
		}

		for(int i=0; i< interfaceElementsEntries.size() ; i++){
			for(int j=0; j < importElements[i].size(); j++){
				if(importElements[i][j] != -1)
					this->edgeElements_->setElementsOfEdgeGlobalEntry(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i]),importElements[i][j]);
			}
		}

	}

}

/*!

 \brief Not all edges are marked with a flag in the beginning. In order to set the correct flags to new points we assign the edge flag of the edge they originated from, similar to the function determineEdgeFlagP2, but this function uses the edgeMap.

*/


template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::assignEdgeFlags(){

	vec2D_LO_Type markedPoints(0);
	EdgeElementsPtr_Type edgeElements = this->getEdgeElements();
	ElementsPtr_Type elements = this->getElementsC();
	MapConstPtr_Type mapRepeatedP1 = this->getMapRepeated();
	vec_int_Type newFlags(edgeElements->numberElements());

	MapConstPtr_Type edgeMap = this->getEdgeMap();

	vec_int_Type markedTrue(edgeElements->numberElements());
	for(int i=0; i< edgeElements->numberElements() ; i++){
		LO p1ID =edgeElements->getElement(i).getNode(0);
		LO p2ID =edgeElements->getElement(i).getNode(1);
		newFlags[i]=this->determineFlagP2(p1ID, p2ID, i , markedPoints );
		if(newFlags[i] != -1){ // questionable point that were given a flag, but that is not certain yet
			vec_LO_Type elementsOfEdge = edgeElements->getElementsOfEdge( (int) i );		
       		for (int j=0; j<elementsOfEdge.size(); j++) {
           		if ( elementsOfEdge[j] == -1 ) 
					markedTrue[i] =1;
			}
		}	
	}
	// It is possible that a Edge is conencted to two surfaces with different Flags, that are also on different Processors
	// This Leads to two different Flags for the same Edge
	// In order to counter that effect we check the interface edges of which we determined the flag via surfaces and check if they have the same flag and if not choose the lower one
	int maxRank = std::get<1>(this->rankRange_);
	if(maxRank >0 ){
		vec_GO_Type edgeSwitch(0);
		vec_LO_Type flags(0);
		for(int i=0; i<edgeElements->numberElements(); i++){
			if(markedTrue[i]==1){
				edgeSwitch.push_back(edgeMap->getGlobalElement(i));
				flags.push_back(newFlags[i]);
			}
		}

		// communticating elements across interface
		Teuchos::ArrayView<GO> edgeSwitchArray = Teuchos::arrayViewFromVector( edgeSwitch);

		MapPtr_Type mapGlobalInterface =
			Teuchos::rcp( new Map_Type( edgeMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgeSwitchArray, 0, this->comm_) );

		// Global IDs of Procs
		// Setting newPoints as to be communicated Values
		MultiVectorLOPtr_Type edgeFlags = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		Teuchos::ArrayRCP< LO > edgeFlagsEntries  = edgeFlags->getDataNonConst(0);

		for(int i=0; i< edgeFlagsEntries.size() ; i++){
			edgeFlagsEntries[i] = flags[i] ;
		}

		MapConstPtr_Type mapGlobalInterfaceUnique = mapGlobalInterface;
		if(mapGlobalInterface->getGlobalNumElements()>0){
			mapGlobalInterfaceUnique = mapGlobalInterface->buildUniqueMap( this->rankRange_ );
		}


		MultiVectorLOPtr_Type isInterfaceElement_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_imp->putScalar( (LO) 0 ); 
		isInterfaceElement_imp->importFromVector( edgeFlags, false, "Insert");

		MultiVectorLOPtr_Type isInterfaceElement_exp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_exp->putScalar( (LO) 0 ); 
		isInterfaceElement_exp->exportFromVector( edgeFlags, false, "Insert");

		MultiVectorLOPtr_Type isInterfaceElement2_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		isInterfaceElement2_imp->putScalar( (LO) 0 ); 
		isInterfaceElement2_imp->importFromVector(isInterfaceElement_imp, false, "Insert");

		isInterfaceElement2_imp->exportFromVector(isInterfaceElement_exp, false, "Insert");

		edgeFlagsEntries  = isInterfaceElement2_imp->getDataNonConst(0);

		for(int i=0; i<edgeFlagsEntries.size(); i++){
			LO entry = edgeMap->getLocalElement(edgeSwitch[i]);
			if(newFlags[entry] > edgeFlagsEntries[i]){
				newFlags[entry] = edgeFlagsEntries[i];

			}
		}
	}
	// In the next Step we need to determine the missing Flags of the 'MarkedMissing' Edges
	// We create a Map of the entries we have and one of the ones we need
	vec_GO_Type edgesNeeded(0); // For the Flag entries we need
	vec_GO_Type edgesActive(0);
	vec_int_Type flagsTmp(0);
	for(int i=0; i<edgeElements->numberElements(); i++){
		if(newFlags[i]==-1){
			edgesNeeded.push_back(edgeMap->getGlobalElement(i));
		}
		else{
			edgesActive.push_back(edgeMap->getGlobalElement(i));
			flagsTmp.push_back(newFlags[i]);
			
		}
	}

	// communticating elements across interface
	Teuchos::ArrayView<GO> edgesNeededArray = Teuchos::arrayViewFromVector( edgesNeeded);
	Teuchos::ArrayView<GO> edgesActiveArray = Teuchos::arrayViewFromVector( edgesActive);
	
	MapPtr_Type mapEdgesNeeded =
		Teuchos::rcp( new Map_Type( edgeMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesNeededArray, 0, this->comm_) );

	MapPtr_Type mapEdgesActive =
		Teuchos::rcp( new Map_Type( edgeMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesActiveArray, 0, this->comm_) );

	MultiVectorLOPtr_Type flagsImport = Teuchos::rcp( new MultiVectorLO_Type( mapEdgesNeeded, 1 ) );
	flagsImport->putScalar(10);

	MultiVectorLOPtr_Type flagsExport = Teuchos::rcp( new MultiVectorLO_Type( mapEdgesActive, 1 ) );
	Teuchos::ArrayRCP< LO > flagExportEntries  = flagsExport->getDataNonConst(0);
	for(int i=0; i< flagExportEntries.size(); i++){
		flagExportEntries[i] = flagsTmp[i];
	}
	
	flagsImport->importFromVector(flagsExport, false, "Insert");

    Teuchos::ArrayRCP< LO > flagImportEntries  = flagsImport->getDataNonConst(0);
	for(int i=0; i<flagImportEntries.size(); i++){
		LO entry = edgeMap->getLocalElement(edgesNeeded[i]);
		if(newFlags[entry] ==-1){
			newFlags[entry] = flagImportEntries[i];
		}
	}

	for(int i=0; i< edgeElements->numberElements() ; i++){
		edgeElements->getElement(i).setFlag(newFlags[i]);
	}
    
}



}
#endif
