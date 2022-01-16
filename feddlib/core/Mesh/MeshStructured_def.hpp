#ifndef MeshStructured_def_hpp
#define MeshStructured_def_hpp
#include "MeshStructured_decl.hpp"

/*!
 Definition of MeshStructured
 
 @brief  MeshStructured
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
//using namespace Teuchos;
namespace FEDD {
template <class SC, class LO, class GO, class NO>
MeshStructured<SC,LO,GO,NO>::MeshStructured():
MeshUnstructured<SC,LO,GO,NO>(),
coorRec(0),
length(0),
height(0),
width(0)
{ }

template <class SC, class LO, class GO, class NO>
MeshStructured<SC,LO,GO,NO>::MeshStructured(CommConstPtrConst_Type& comm):
MeshUnstructured<SC,LO,GO,NO>(comm),
coorRec(0),
length(0),
height(0),
width(0)
{ }

template <class SC, class LO, class GO, class NO>
MeshStructured<SC,LO,GO,NO>::~MeshStructured(){

}

template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::setGeometry2DRectangle(std::vector<double> coordinates, double l, double h){
    
	coorRec	= coordinates;
	length 	= l;
	height 	= h;
    
    this->dim_ = 2;
}

template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::setGeometry3DBox(std::vector<double> coordinates, double l, double w, double h){
    
    coorRec	= coordinates;
    length 	= l;
    width 	= w;
    height 	= h;
    
    this->dim_ = 3;
}
    
template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::setRankRange(int numProcsCoarseSolve){
    get<0>(this->rankRange_) = 0;
    get<1>(this->rankRange_) = this->comm_->getSize() - 1 - numProcsCoarseSolve;
}

template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::buildMesh2DTPM(std::string FEType,
                                              int N,
                                              int M,
                                              int numProcsCoarseSolve,
                                              std::string underlyingLib){
    
    buildMesh2D( FEType, N, M, numProcsCoarseSolve, underlyingLib );
    
    setRankRange( numProcsCoarseSolve );
    
    buildSurfaceLinesSquare();

	// New Feature that adds edgeFlags and surfaceFlags
	// setEdgeFlags
	// setSurfaceFlags
        
}
    
template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::buildMesh2DMiniTPM(std::string FEType,
                                                     int N,
                                                     int M,
                                                     int numProcsCoarseSolve,
                                                     std::string underlyingLib){
    
    this->FEType_ = FEType;
        
    this->numElementsGlob_ = 4;
    int nmbPoints;

    vec2D_int_ptr_Type elementsVec;
    if (FEType=="P2") {
        nmbPoints = 15;
        elementsVec.reset(new std::vector<std::vector<int> >(this->numElementsGlob_,std::vector<int>(6,-1)));
    } else if(FEType=="P1"){
        nmbPoints = 6;
        elementsVec.reset(new std::vector<std::vector<int> >(this->numElementsGlob_,std::vector<int>(3,-1)));
    }

    this->pointsRep_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(2,0.0)));
    this->bcFlagRep_.reset(new std::vector<int> (nmbPoints,0));
    this->pointsUni_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(2,0.0)));
    this->bcFlagUni_.reset(new std::vector<int> (nmbPoints,0));
    
    Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);
    for (int i=0; i<nmbPoints; i++) {
        pointsRepGlobMapping[i] = i;
    }
    
    this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );

    this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );
    
    double h = 0.1;
    int counter = 0;
    if (FEType=="P2") {
        for (int i=0; i<3; i++) {
            for (int j=0; j<5; j++) {
                (*this->pointsRep_)[counter][0] = j*h/2;
                (*this->pointsRep_)[counter][1] = i*h/2;
                
                (*this->pointsUni_)[counter][0] = j*h/2;
                (*this->pointsUni_)[counter][1] = i*h/2;
                counter++;
            }
        }
    } else if(FEType=="P1"){
        for (int i=0; i<2; i++) {
            for (int j=0; j<3; j++) {
                (*this->pointsRep_)[counter][0] = j*h;
                (*this->pointsRep_)[counter][1] = i*h;
                
                (*this->pointsUni_)[counter][0] = j*h;
                (*this->pointsUni_)[counter][1] = i*h;
                counter++;
            }
        }
    }

    vec_int_ptr_Type elementFlag = Teuchos::rcp( new vec_int_Type( elementsVec->size(), 0 ) );
    
    counter = 0;
    int S=1;
    int R=2;
    int P2M = 2*(R+1)-1;
    if (FEType=="P2") {
        
         for (int s=0; s < S; s++) {
                for (int r=0; r < R; r++) {
                        
                    (*elementsVec)[counter][0] = 2*(r+1)    + 2*P2M * (s) ;
                    (*elementsVec)[counter][1] = 2*(r)      + 2*P2M * (s) ;
                    (*elementsVec)[counter][2] = 2*(r+1)    + 2*P2M * (s+1) ;
                    
                    (*elementsVec)[counter][3] = 2*(r) +1    + 2*P2M * (s) ;
                    (*elementsVec)[counter][4] = 2*(r) +1    + 2*P2M * (s) +P2M ;
                    (*elementsVec)[counter][5] = 2*(r+1)    + 2*P2M * (s) +P2M ;

                    counter++;
                    
                    
                    
                    (*elementsVec)[counter][0] = 2*(r)     + 2*P2M * (s+1) ;
                    (*elementsVec)[counter][1] = 2*(r)     + 2*P2M * (s) ;
                    (*elementsVec)[counter][2] = 2*(r+1)    + 2*P2M * (s+1) ;
                    
                    (*elementsVec)[counter][3] = 2*(r)        + 2*P2M * (s) +P2M ;
                    (*elementsVec)[counter][4] = 2*(r) +1     + 2*P2M * (s) +P2M ;
                    (*elementsVec)[counter][5] = 2*(r) +1    + 2*P2M * (s+1) ;
                    
                    counter++;
                }
         }
        
    } else if(FEType=="P1") {
       
        for (int s=0; s < S; s++) {
            for (int r=0; r < R; r++) {
                
                (*elementsVec)[counter][0] = r+1 + (R+1)* s;
                (*elementsVec)[counter][1] = r + (R+1)* s;
                (*elementsVec)[counter][2] = r+1 + (R+1) * (s+1);
                
                counter++;
                
                (*elementsVec)[counter][0] = r + (R+1) * (s+1);
                (*elementsVec)[counter][1] = r + (R+1) * (s);
                (*elementsVec)[counter][2] = r+1 + (R+1) * (s+1);

                counter++;
            }
            
        }
    }
    
    setRankRange( numProcsCoarseSolve );
    
    buildElementsClass( elementsVec, elementFlag  );
    
    buildSurfaceLinesSquareMiniTPM( FEType );
    
}


template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::buildElementsClass( vec2D_int_ptr_Type elements, vec_int_ptr_Type elementFlag ){

    this->elementsC_.reset(new Elements ( this->FEType_, this->dim_ ) );
    bool setFlags = !elementFlag.is_null();
    for (int i=0; i<elements->size(); i++) {
        std::vector<LO> tmpElement;
        for (int j=0; j<elements->at(i).size(); j++) {
            tmpElement.push_back( (*elements)[i][j] );
        }
        if (setFlags){
            FiniteElement fe( tmpElement, (*elementFlag)[i] );
            this->elementsC_->addElement( fe );
        }
        else{
            FiniteElement fe( tmpElement );
            this->elementsC_->addElement( fe );
        }
    }
}
    
template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::buildSurfaceLinesSquare(){

//    for (int i=0; i<this->elementsC_->numberElements(); i++) {
//
//    }
    
}

template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::buildSurfaceLinesSquareMiniTPM(string feType){
    ElementsPtr_Type elementsMesh = this->getElementsC();
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "Must be implemented for new elements!");
    if (feType=="P2") {
//        vec_int_Type tmpSurface(3);
//        tmpSurface[0] = 0; tmpSurface[1] = 2; tmpSurface[2] = 1;
//        elementsMesh->getElement(0).setLocalSurface( 0, tmpSurface, 1 );
//
//        tmpSurface[0] = 0; tmpSurface[1] = 10; tmpSurface[2] = 5;
//        elementsMesh->getElement(1).setLocalSurface( 1, tmpSurface, 2 );
//        tmpSurface[0] = 10; tmpSurface[1] = 12; tmpSurface[2] = 11;
//        elementsMesh->getElement(1).setLocalSurface( 2, tmpSurface, 1 );
//
//
//        tmpSurface[0] = 2; tmpSurface[1] = 4; tmpSurface[2] = 3;
//        elementsMesh->getElement(2).setLocalSurface( 3, tmpSurface, 1 );
//        tmpSurface[0] = 4; tmpSurface[1] = 14; tmpSurface[2] = 9;
//        elementsMesh->getElement(2).setLocalSurface( 4, tmpSurface, 3 );
//
//        tmpSurface[0] = 12; tmpSurface[1] = 14; tmpSurface[2] = 13;
//        elementsMesh->getElement(3).setLocalSurface( 5, tmpSurface, 1 );
        
    } else {
//        vec_int_Type tmpSurface(2);
//        tmpSurface[0] = 0; tmpSurface[1] = 1;
//        elementsMesh->getElement(0).setLocalSurface( 0, tmpSurface, 1 );
//
//        tmpSurface[0] = 0; tmpSurface[1] = 3;
//        elementsMesh->getElement(1).setLocalSurface( 1, tmpSurface, 2 );
//        tmpSurface[0] = 3; tmpSurface[1] = 4;
//        elementsMesh->getElement(1).setLocalSurface( 2, tmpSurface, 1 );
//
//        tmpSurface[0] = 1; tmpSurface[1] = 2;
//        elementsMesh->getElement(2).setLocalSurface( 3, tmpSurface, 1 );
//        tmpSurface[0] = 2; tmpSurface[1] = 5;
//        elementsMesh->getElement(2).setLocalSurface( 4, tmpSurface, 3 );
//
//        tmpSurface[0] = 4; tmpSurface[1] = 5;
//        elementsMesh->getElement(3).setLocalSurface( 5, tmpSurface, 1 );
    }
    
}

template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::buildMesh2D(std::string FEType,
                                                 int N,
                                                 int M,
                                                 int numProcsCoarseSolve,
                                                 std::string underlyingLib){
    
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ScalarTraits;
    
    
    TEUCHOS_TEST_FOR_EXCEPTION(!(M>=1),std::logic_error,"H/h is to small.");
    TEUCHOS_TEST_FOR_EXCEPTION(this->comm_.is_null(),std::runtime_error,"comm_ is null.");
    
    bool verbose (this->comm_->getRank() == 0);
    
    setRankRange( numProcsCoarseSolve );
    
    if (verbose) {
        cout << endl;
    }
    
    int         rank = this->comm_->getRank();
    int         size = this->comm_->getSize() - numProcsCoarseSolve;
    
    SC      h = length/(M*N);//add variable length/width/heigth
    SC      H = length/N;
    GO 	nmbPoints_oneDir;
    
    LO nmbElements;
    LO nmbPoints;
    vec2D_int_ptr_Type elementsVec;
    vec_int_ptr_Type elementFlag;
    
    if (FEType == "P0") {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"implement P0.");
    }
    else if (FEType == "P1") {
        nmbPoints_oneDir 	= N * (M+1) - (N-1);
        nmbPoints			= (M+1)*(M+1);
    }
    else if(FEType == "P2"){
        nmbPoints_oneDir 	=  N * (2*(M+1)-1) - (N-1);
        nmbPoints 			= (2*(M+1)-1)*(2*(M+1)-1);
    }
    else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Wrong FE-Type, either P1 or P2.");
    }
    
    this->FEType_ = FEType;
    
    GO nmbPGlob_oneDir = N * (M+1) - (N-1);
    
    this->numElementsGlob_ = 2*(nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1);

    if (rank>=size) {
        M = -1; // keine Schleife wird ausgefuehrt
        nmbElements = 0;
        nmbPoints = 0;
    }
    else{
        nmbElements = 2*(M)*(M);
    }
        
    // P1 Mesh
    if (FEType == "P1") {
        if (verbose) {
            cout << "-- H:"<<H << " h:" <<h << " --" << endl;
        }
        if (verbose) {
            cout << "-- Building P1 Points Repeated ... " << endl;
        }
        
        this->pointsRep_.reset(new vec2D_dbl_Type(nmbPoints,std::vector<double>(2,0.0)));
        this->bcFlagRep_.reset(new vec_int_Type (nmbPoints,0));
        elementsVec = Teuchos::rcp(new vec2D_int_Type(nmbElements,std::vector<int>(3,-1)));
        Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);
        
        int counter = 0;
        int offset_x = (rank % N);
        int offset_y = 0;
        
        if ((rank % (N*N))>=N) {
            offset_y = (int) (rank % (N*N))/(N);
        }
        
        for (int s=0; s < M+1; s++) {
            for (int r=0; r < M+1; r++) {
                (*this->pointsRep_)[counter][0] = r*h + offset_x * H;
                if ((*this->pointsRep_)[counter][0]<100*ScalarTraits<SC>::eps() && (*this->pointsRep_)[counter][0]>-100*ScalarTraits<SC>::eps()) { (*this->pointsRep_)[counter][0]=0.0;}
                (*this->pointsRep_)[counter][1] = s*h + offset_y * H;
                if ((*this->pointsRep_)[counter][1]<100*ScalarTraits<SC>::eps() && (*this->pointsRep_)[counter][1]>-100*ScalarTraits<SC>::eps()) {(*this->pointsRep_)[counter][1]=0.0;}
                pointsRepGlobMapping[counter] = r + s*nmbPoints_oneDir + offset_x*(M) + offset_y*(nmbPoints_oneDir)*M;
                if ((*this->pointsRep_)[counter][0] > (coorRec[0]+length-100*ScalarTraits<SC>::eps()) 	|| (*this->pointsRep_)[counter][0] < (coorRec[0]+100*ScalarTraits<SC>::eps()) ||
                    (*this->pointsRep_)[counter][1] > (coorRec[1]+height-100*ScalarTraits<SC>::eps()) 	|| (*this->pointsRep_)[counter][1] < (coorRec[1]+100*ScalarTraits<SC>::eps())) {
                    
                    (*this->bcFlagRep_)[counter] = 1;
                }
                counter++;
            }
        }
        
        if (verbose) {
            cout << " done! --" << endl;
        }
        
        if (verbose) {
            cout << "-- Building P1 Repeated and Unique Map ... " << flush;
        }

        this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );

        this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );

        if (verbose) {
            cout << " done! --" << endl;
        }
        
        if (verbose) {
            cout << "-- Building P1 Unique Points ... " << flush;
        }
        
        this->pointsUni_.reset(new std::vector<std::vector<double> >(this->mapUnique_->getNodeNumElements(),std::vector<double>(2,0.0)));
        this->bcFlagUni_.reset(new std::vector<int> (this->mapUnique_->getNodeNumElements(),0));
        LO index;
        for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
            
            index = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement(i) );
            
            (*this->pointsUni_)[i][0] = (*this->pointsRep_)[index][0];
            (*this->pointsUni_)[i][1] = (*this->pointsRep_)[index][1];
            (*this->bcFlagUni_)[i] = (*this->bcFlagRep_)[index];
        }
        if (verbose) {
            cout << " done! --" << endl;
        }
        
        
        if (verbose) {
            cout << "-- Building P1 Elements ... " << flush;
        }
        vec_int_ptr_Type elementFlag = Teuchos::rcp(new vec_int_Type( elementsVec->size(),0 ) );
        counter = 0;
        double x_ref, y_ref;
        
        for (int s=0; s < M; s++) {
            for (int r=0; r < M; r++) {
                
                (*elementsVec)[counter][0] = r+1 + (M+1) * s;
                (*elementsVec)[counter][1] = r + (M+1) * s ;
                (*elementsVec)[counter][2] = r+1 + (M+1) * (s+1);
                x_ref = ( this->pointsRep_->at( elementsVec->at(counter).at(0) ).at(0) + this->pointsRep_->at( elementsVec->at(counter).at(1) ).at(0) + this->pointsRep_->at( elementsVec->at(counter).at(2) ).at(0) ) / 3.;
                y_ref = ( this->pointsRep_->at( elementsVec->at(counter).at(0) ).at(1) + this->pointsRep_->at( elementsVec->at(counter).at(1) ).at(1) + this->pointsRep_->at( elementsVec->at(counter).at(2) ).at(1) ) / 3.;
                if ( x_ref>=0.3  && x_ref<=0.7) {
                    if ( y_ref>= 0.6) {
                        elementFlag->at(counter) = 1;
                    }
                }

                counter++;
                
                (*elementsVec)[counter][0] = r + (M+1) * (s+1);
                (*elementsVec)[counter][1] = r + (M+1) * (s);
                (*elementsVec)[counter][2] = r+1 + (M+1) * (s+1);

                x_ref = ( this->pointsRep_->at( elementsVec->at(counter).at(0) ).at(0) + this->pointsRep_->at( elementsVec->at(counter).at(1) ).at(0) + this->pointsRep_->at( elementsVec->at(counter).at(2) ).at(0) ) / 3.;
                y_ref = ( this->pointsRep_->at( elementsVec->at(counter).at(0) ).at(1) + this->pointsRep_->at( elementsVec->at(counter).at(1) ).at(1) + this->pointsRep_->at( elementsVec->at(counter).at(2) ).at(1) ) / 3.;
                if ( x_ref>=0.3  && x_ref<=0.7) {
                    if ( y_ref>= 0.6) {
                        elementFlag->at(counter) = 1;
                    }
                }

                counter++;
            }
        }
        
        if (verbose) {
            cout << " done! --" << endl;
        }
    }
    
    // P2 Mesh
    
    
    else if(FEType == "P2"){
        if (verbose) {
            cout << "-- H:"<<H << " h:" <<h << " --" << endl;
        }
        if (verbose) {
            cout << "-- Building P2 Points Repeated ... " << flush;
        }
        
        this->pointsRep_.reset(new vec2D_dbl_Type(nmbPoints, vec_dbl_Type(2,0.0)));
        this->bcFlagRep_.reset(new vec_int_Type (nmbPoints,0));
        elementsVec = Teuchos::rcp(new vec2D_int_Type(nmbElements, vec_int_Type(6,-1)));
        Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);
        
        int counter = 0;
        int offset_x = (rank % N);
        int offset_y = 0;
        
        if ((rank % (N*N))>=N) {
            offset_y = (int) (rank % (N*N))/(N);
        }
        bool p1point;
        int p1_s = 0;
        int p1_r = 0;
        for (int s=0; s < 2*(M+1)-1; s++) {
            for (int r=0; r < 2*(M+1)-1; r++) {
                p1point = false;
                if (s%2==0 && r%2==0) {
                    p1point = true;
                    p1_s = s/2;
                    p1_r = r/2;
                }
                (*this->pointsRep_)[counter][0] = r*h/2.0 + offset_x * H;
                if ((*this->pointsRep_)[counter][0]<100*ScalarTraits<SC>::eps() && (*this->pointsRep_)[counter][0]>-100*ScalarTraits<SC>::eps()) (*this->pointsRep_)[counter][0]=0.0;
                (*this->pointsRep_)[counter][1] = s*h/2.0 + offset_y * H;
                if ((*this->pointsRep_)[counter][1]<100*ScalarTraits<SC>::eps() && (*this->pointsRep_)[counter][1]>-100*ScalarTraits<SC>::eps()) (*this->pointsRep_)[counter][1]=0.0;
                
                pointsRepGlobMapping[counter] = r + s*nmbPoints_oneDir + offset_x*(2*(M+1)-2) + offset_y*(nmbPoints_oneDir)*(2*(M+1)-2) ;
                
                if ((*this->pointsRep_)[counter][0] > (coorRec[0]+length-100*ScalarTraits<SC>::eps()) 	|| (*this->pointsRep_)[counter][0] < (coorRec[0]+100*ScalarTraits<SC>::eps()) ||
                    (*this->pointsRep_)[counter][1] > (coorRec[1]+height-100*ScalarTraits<SC>::eps()) 	|| (*this->pointsRep_)[counter][1] < (coorRec[1]+100*ScalarTraits<SC>::eps())){
                    (*this->bcFlagRep_)[counter] = 1;
                    
                }
                counter++;
            }
        }
        
        if (verbose)
            cout << " done! --" << endl;
        
        if (verbose)
            cout << "-- Building P2 Repeated and Unique Map ... " << flush;
        
        this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );
        
        this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );
        
        if (verbose) {
            cout << " done! --" << endl;
        }
        
        if (verbose) {
            cout << "-- Building P2 Unique Points ... " << flush;
        }
        
        this->pointsUni_.reset(new std::vector<std::vector<double> >(this->mapUnique_->getNodeNumElements(),std::vector<double>(2,0.0)));
        this->bcFlagUni_.reset(new std::vector<int> (this->mapUnique_->getNodeNumElements(),0));
        
        LO index;
        for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
            
            index = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement(i) );
            
            (*this->pointsUni_)[i][0] = (*this->pointsRep_)[index][0];
            (*this->pointsUni_)[i][1] = (*this->pointsRep_)[index][1];
            (*this->bcFlagUni_)[i] = (*this->bcFlagRep_)[index];
        }
        
        if (verbose) {
            cout << " done! --" << endl;
        }
        
        //                Triangle numbering
        //                    2
        //                  * *
        //                *   *
        //              4	  5
        //            *       *
        //          *         *
        //        1 * * 3 * * 0
        
        
        if (verbose)
            cout << "-- Building P2 Elements ... " << flush;
        
        int    P2M = 2*(M+1)-1;
        
        vec_int_ptr_Type elementFlag = Teuchos::rcp(new vec_int_Type( elementsVec->size(),0 ) );
        counter = 0;
        double x_ref, y_ref;
        
        for (int s=0; s < M; s++) {
            for (int r=0; r < M; r++) {
                
                (*elementsVec)[counter][0] = 2*(r+1)	+ 2*P2M * (s) ;
                (*elementsVec)[counter][1] = 2*(r)      + 2*P2M * (s) ;
                (*elementsVec)[counter][2] = 2*(r+1)	+ 2*P2M * (s+1) ;
                
                (*elementsVec)[counter][3] = 2*(r) +1	+ 2*P2M * (s) ;
                (*elementsVec)[counter][4] = 2*(r) +1	+ 2*P2M * (s) +P2M ;
                (*elementsVec)[counter][5] = 2*(r+1)	+ 2*P2M * (s) +P2M ;

                x_ref = ( this->pointsRep_->at( elementsVec->at(counter).at(0) ).at(0) + this->pointsRep_->at( elementsVec->at(counter).at(1) ).at(0) + this->pointsRep_->at( elementsVec->at(counter).at(2) ).at(0) ) / 3.;
                y_ref = ( this->pointsRep_->at( elementsVec->at(counter).at(0) ).at(1) + this->pointsRep_->at( elementsVec->at(counter).at(1) ).at(1) + this->pointsRep_->at( elementsVec->at(counter).at(2) ).at(1) ) / 3.;
                if ( x_ref>=0.3  && x_ref<=0.7) {
                    if ( y_ref>= 0.6) {
                        elementFlag->at(counter) = 1;
                    }
                }
                
                counter++;
                
                (*elementsVec)[counter][0] = 2*(r) 	+ 2*P2M * (s+1) ;
                (*elementsVec)[counter][1] = 2*(r) 	+ 2*P2M * (s) ;
                (*elementsVec)[counter][2] = 2*(r+1)	+ 2*P2M * (s+1) ;
                
                (*elementsVec)[counter][3] = 2*(r)		+ 2*P2M * (s) +P2M ;
                (*elementsVec)[counter][4] = 2*(r) +1 	+ 2*P2M * (s) +P2M ;
                (*elementsVec)[counter][5] = 2*(r) +1	+ 2*P2M * (s+1) ;
                
                x_ref = ( this->pointsRep_->at( elementsVec->at(counter).at(0) ).at(0) + this->pointsRep_->at( elementsVec->at(counter).at(1) ).at(0) + this->pointsRep_->at( elementsVec->at(counter).at(2) ).at(0) ) / 3.;
                y_ref = ( this->pointsRep_->at( elementsVec->at(counter).at(0) ).at(1) + this->pointsRep_->at( elementsVec->at(counter).at(1) ).at(1) + this->pointsRep_->at( elementsVec->at(counter).at(2) ).at(1) ) / 3.;
                if ( x_ref>=0.3  && x_ref<=0.7) {
                    if ( y_ref>= 0.6) {
                        elementFlag->at(counter) = 1;
                    }
                }
                
                counter++;
                
            }
        }
        
        
        
        if (verbose) {
            cout << " done! --" << endl;
        }
    }
	// At the end of these functions we have created the node and element list and set the corresponding flags. As mesh structured is derived from the mesh class, edges and surface are not set yet. 
	// Consequentialy there are no subelements set for elements and thus no surfaceIntegrals can be assembled

    buildElementsClass(elementsVec, elementFlag);
 
}
    
template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::buildMesh3D(std::string FEType,
                                                 int N,
                                                 int M,
                                                 int numProcsCoarseSolve,
                                                 std::string underlyingLib){
    
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ScalarTraits;
    
    TEUCHOS_TEST_FOR_EXCEPTION(!(M>=1),std::logic_error,"H/h is to small.");
    TEUCHOS_TEST_FOR_EXCEPTION(this->comm_.is_null(),std::runtime_error,"comm_ is null.");
    
    bool verbose (this->comm_->getRank() == 0);
    
    setRankRange( numProcsCoarseSolve );
    
    if (verbose) {
        cout << endl;
    }

    SC eps = ScalarTraits<SC>::eps();
    
    int         rank = this->comm_->getRank();
    int         size = this->comm_->getSize() - numProcsCoarseSolve;
    
    SC      h = length/(M*N);//add variable length/width/heigth
    SC      H = length/N;
   
    LO 	nmbPoints_oneDir;
    
    LO nmbElements;
    LO nmbPoints;
    if (FEType == "P0") {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"implement P0.");
    }
    else if (FEType == "P1") {
        nmbPoints_oneDir 	= N * (M+1) - (N-1);
        nmbPoints			= (M+1)*(M+1)*(M+1);
    }
    else if(FEType == "P2"){
        nmbPoints_oneDir 	=  N * (2*(M+1)-1) - (N-1);
        nmbPoints 			= (2*(M+1)-1)*(2*(M+1)-1)*(2*(M+1)-1);
    }
    else if(FEType == "P2-CR"){
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"P2-CR might not work properly.");
        nmbPoints_oneDir 	=  N * (2*(M+1)-1) - (N-1);
        nmbPoints 			= (2*(M+1)-1)*(2*(M+1)-1)*(2*(M+1)-1);
    }
    else if(FEType == "P1-disc" || FEType == "P1-disc-global"){
        
    }
    else if(FEType == "Q1"){

    }
    else if(FEType == "Q2"){
     
    }
    else if(FEType == "Q2-20"){
        
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Wrong FE-Type, only P1,P1-disc, P1-disc-global, P2, P2-CR, Q1, Q2, Q2-20.");

    this->FEType_ = FEType;
    
    GO nmbPGlob_oneDir = N * (M+1) - (N-1);
    
    this->numElementsGlob_ = 6*(nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1);
    int MM=M;
    if (rank>=size) {
        M = -1; // keine Schleife wird ausgefuehrt
        nmbElements = 0;
        nmbPoints = 0;
    }
    else{
        nmbElements = 6*M*M*M;
    }
    vec2D_int_ptr_Type elementsVec;
    vec_int_ptr_Type elementFlag;
    // P1 Mesh
    if (FEType == "P1") {
        
        this->pointsRep_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(3,0.0)));
        this->bcFlagRep_.reset(new std::vector<int> (nmbPoints,0));
        elementsVec = Teuchos::rcp(new vec2D_int_Type( nmbElements, vec_int_Type(4, -1) ));
        elementFlag = Teuchos::rcp(new vec_int_Type( elementsVec->size(),0 ) );
        Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);
        
        int counter = 0;
        int offset_x = (rank % N);
        int offset_y = 0;
        int offset_z = 0;
        
        if ((rank % (N*N))>=N) {
            offset_y = (int) (rank % (N*N))/(N);
        }
        
        if ((rank % (N*N*N))>=N*N ) {
            offset_z = (int) (rank % (N*N*N))/(N*(N));
        }
        
        for (int t=0; t < M+1; t++) {
            for (int s=0; s < M+1; s++) {
                for (int r=0; r < M+1; r++) {
                    (*this->pointsRep_)[counter][0] = r*h + offset_x * H;
                    if ((*this->pointsRep_)[counter][0]<eps && (*this->pointsRep_)[counter][0]>-eps) (*this->pointsRep_)[counter][0]=0.0;
                    
                    (*this->pointsRep_)[counter][1] = s*h + offset_y * H;
                    if ((*this->pointsRep_)[counter][1]<eps && (*this->pointsRep_)[counter][1]>-eps) (*this->pointsRep_)[counter][1]=0.0;
                    
                    (*this->pointsRep_)[counter][2] = t*h + offset_z * H;
                    if ((*this->pointsRep_)[counter][2]<eps && (*this->pointsRep_)[counter][2]>-eps) (*this->pointsRep_)[counter][2]=0.0;
                    
                    pointsRepGlobMapping[counter] = r + s*nmbPoints_oneDir + t*nmbPoints_oneDir*nmbPoints_oneDir \
                    + offset_x*(M) + offset_y*(nmbPoints_oneDir)*M + offset_z*(nmbPoints_oneDir)*(nmbPoints_oneDir)*M  ;
                    
                    if ((*this->pointsRep_)[counter][0] > (coorRec[0]+length-eps) 	|| (*this->pointsRep_)[counter][0] < (coorRec[0]+eps) ||
                        (*this->pointsRep_)[counter][1] > (coorRec[1]+width-eps) 	|| (*this->pointsRep_)[counter][1] < (coorRec[1]+eps) ||
                        (*this->pointsRep_)[counter][2] > (coorRec[2]+height-eps) 	|| (*this->pointsRep_)[counter][2] < (coorRec[2]+eps) ) {
                        
                        (*this->bcFlagRep_)[counter] = 1;
                    }
                    
                    counter++;
                }
            }
        }
        
        this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );

        this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );
        
        this->pointsUni_.reset(new std::vector<std::vector<double> >(this->mapUnique_->getNodeNumElements(),std::vector<double>(3,0.0)));
        this->bcFlagUni_.reset(new std::vector<int> (this->mapUnique_->getNodeNumElements(),0));
        LO index;
        for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
            
            index = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement(i) );
            
            (*this->pointsUni_)[i][0] = (*this->pointsRep_)[index][0];
            (*this->pointsUni_)[i][1] = (*this->pointsRep_)[index][1];
            (*this->pointsUni_)[i][2] = (*this->pointsRep_)[index][2];
            (*this->bcFlagUni_)[i] = (*this->bcFlagRep_)[index];
        }
        
        counter = 0;
        for (int t=0; t < M; t++) {
            for (int s=0; s < M; s++) {
                for (int r=0; r < M; r++) {
                    (*elementsVec)[counter][0] = r+1 + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][1] = r + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][2] = r+1 + (M+1) * s + (M+1)*(M+1) * (t+1) ;
                    (*elementsVec)[counter][3] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    counter++;
                    (*elementsVec)[counter][0] = r + (M+1) * s + (M+1)*(M+1) * (t+1) ;
                    (*elementsVec)[counter][1] = r + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][2] = r+1 + (M+1) * s + (M+1)*(M+1) * (t+1) ;
                    (*elementsVec)[counter][3] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    counter++;
                    (*elementsVec)[counter][0] = r+1 + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][1] = r + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][2] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][3] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    counter++;
                    (*elementsVec)[counter][0] = r + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][1] = r + (M+1) * (s+1) + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][2] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][3] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    counter++;
                    (*elementsVec)[counter][0] = r + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][1] = r + (M+1) * (s+1) + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][2] = r + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    (*elementsVec)[counter][3] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    counter++;
                    (*elementsVec)[counter][0] = r + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][1] = r + (M+1) * s + (M+1)*(M+1) * (t+1) ;
                    (*elementsVec)[counter][2] = r + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    (*elementsVec)[counter][3] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    counter++;
                }
            }
        }
        buildElementsClass(elementsVec, elementFlag);
    }
    
    else if(FEType == "P2"){
        
        this->pointsRep_.reset(new vec2D_dbl_Type(nmbPoints, vec_dbl_Type(3, 0.0)));
        this->bcFlagRep_.reset(new vec_int_Type (nmbPoints, 0));
        elementsVec = Teuchos::rcp(new vec2D_int_Type(nmbElements, vec_int_Type(10, -1)));
        elementFlag = Teuchos::rcp(new vec_int_Type( elementsVec->size(),0 ) );
        Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);
        
        int counter = 0;
        int offset_x = (rank % N);
        int offset_y = 0;
        int offset_z = 0;
        
        if ((rank % (N*N))>=N) {
            offset_y = (int) (rank % (N*N))/(N);
        }
        
        if ((rank % (N*N*N))>=N*N ) {
            offset_z = (int) (rank % (N*N*N))/(N*(N));
        }
        bool p1point;
        int p1_s = 0;
        int p1_r = 0;
        int p1_t = 0;
        for (int t=0; t < 2*(M+1)-1; t++) {
            for (int s=0; s < 2*(M+1)-1; s++) {
                for (int r=0; r < 2*(M+1)-1; r++) {
                    p1point = false;
                    if (s%2==0 && r%2==0 && t%2==0) {
                        p1point = true;
                        p1_s = s/2;
                        p1_r = r/2;
                        p1_t = t/2;
                    }
                    (*this->pointsRep_)[counter][0] = r*h/2.0 + offset_x * H;
                    if ((*this->pointsRep_)[counter][0]<eps && (*this->pointsRep_)[counter][0]>-eps) (*this->pointsRep_)[counter][0]=0.0;
                    (*this->pointsRep_)[counter][1] = s*h/2.0 + offset_y * H;
                    if ((*this->pointsRep_)[counter][1]<eps && (*this->pointsRep_)[counter][1]>-eps) (*this->pointsRep_)[counter][1]=0.0;
                    (*this->pointsRep_)[counter][2] = t*h/2.0 + offset_z * H;
                    if ((*this->pointsRep_)[counter][2]<eps && (*this->pointsRep_)[counter][2]>-eps) (*this->pointsRep_)[counter][2]=0.0;
                    
                    pointsRepGlobMapping[counter] = r + s*nmbPoints_oneDir + t*nmbPoints_oneDir*nmbPoints_oneDir \
                    + offset_x*(2*(M+1)-2) + offset_y*(nmbPoints_oneDir)*(2*(M+1)-2) + offset_z*(nmbPoints_oneDir)*(nmbPoints_oneDir)*(2*(M+1)-2) ;
                    
                    if ((*this->pointsRep_)[counter][0] > (coorRec[0]+length-eps) || (*this->pointsRep_)[counter][0] < (coorRec[0]+eps) ||
                        (*this->pointsRep_)[counter][1] > (coorRec[1]+width-eps) 	|| (*this->pointsRep_)[counter][1] < (coorRec[1]+eps) ||
                        (*this->pointsRep_)[counter][2] > (coorRec[2]+height-eps) || (*this->pointsRep_)[counter][2] < (coorRec[2]+eps) ) {
                        (*this->bcFlagRep_)[counter] = 1;
                        
                    }
                    
                    counter++;
                }
            }
        }
        
        this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );
        
        this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );
        
        this->pointsUni_.reset(new std::vector<std::vector<double> >(this->mapUnique_->getNodeNumElements(),std::vector<double>(3,0.0)));
        this->bcFlagUni_.reset(new std::vector<int> (this->mapUnique_->getNodeNumElements(),0));
        
        for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
            (*this->pointsUni_)[i][0] = (this->mapUnique_->getGlobalElement(i) % nmbPoints_oneDir) * h/2;
            if ((*this->pointsUni_)[i][0]<eps && (*this->pointsUni_)[i][0]>-eps) (*this->pointsUni_)[i][0]=0.0;
            
            (*this->pointsUni_)[i][1] = ((int) ((this->mapUnique_->getGlobalElement(i) % (nmbPoints_oneDir*nmbPoints_oneDir)) / nmbPoints_oneDir) + eps) *h/2;
            if ((*this->pointsUni_)[i][1]<eps && (*this->pointsUni_)[i][1]>-eps) (*this->pointsUni_)[i][1]=0.0;
            
            (*this->pointsUni_)[i][2] = ((int)(this->mapUnique_->getGlobalElement(i) / (nmbPoints_oneDir*nmbPoints_oneDir) + eps)) * h/2;
            if ((*this->pointsUni_)[i][2]<eps && (*this->pointsUni_)[i][2]>-eps) (*this->pointsUni_)[i][2]=0.0;
            
            if ((*this->pointsUni_)[i][0] > (coorRec[0]+length-eps) 	|| (*this->pointsUni_)[i][0] < (coorRec[0]+eps) ||
                (*this->pointsUni_)[i][1] > (coorRec[1]+width-eps) 	|| (*this->pointsUni_)[i][1] < (coorRec[1]+eps) ||
                (*this->pointsUni_)[i][2] > (coorRec[2]+height-eps) 	|| (*this->pointsUni_)[i][2] < (coorRec[2]+eps) ) {
                (*this->bcFlagUni_)[i] = 1;
                
            }
        }
        
        //                Face 1          Face2               Face 3            Face 4
        //                    2      2 * * 9 * * 3        3 * * 9 * * 2          	3
        //                  * *      *          *          *          * 		  * *
        //                *   *      *        *             *        *          *   *
        //              5	  6      6      7                8      5         8	    7
        //            *       *      *    *                   *    *        *       *
        //          *         *      *  *                      *  *       *         *
        //        1 * * 4 * * 0       0                         1       1 * * 4 * * 0
        
        
        int    P2M = 2*(M+1)-1;
        
        counter = 0;
        for (int t=0; t < M; t++) {
            for (int s=0; s < M; s++) {
                for (int r=0; r < M; r++) {
                    
                    (*elementsVec)[counter][0] = 2*(r+1)	+ 2*P2M * (s) 	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][1] = 2*(r) 	+ 2*P2M * (s) 	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][2] = 2*(r+1)	+ 2*P2M * (s) 	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][3] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    
                    (*elementsVec)[counter][4] = 2*(r) +1	+ 2*P2M * (s)		+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][6] = 2*(r+1)		+ 2*P2M * (s)		+ 2*P2M*P2M * (t) +P2M*P2M;
                    (*elementsVec)[counter][7] = 2*(r+1)		+ 2*P2M * (s) +P2M 	+ 2*P2M*P2M * (t) +P2M*P2M;
                    (*elementsVec)[counter][5] = 2*(r) +1	+ 2*P2M * (s) 		+ 2*P2M*P2M * (t) +P2M*P2M;
                    (*elementsVec)[counter][8] = 2*(r) +1	+ 2*P2M * (s) +P2M 	+ 2*P2M*P2M * (t) +P2M*P2M;
                    (*elementsVec)[counter][9] = 2*(r+1)		+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t+1);
                    
                    counter++;
                    
                    (*elementsVec)[counter][0] = 2*(r) 	+ 2*P2M * (s)	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][1] = 2*(r) 	+ 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][2] = 2*(r+1)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][3] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    
                    (*elementsVec)[counter][4] = 2*(r) 		+ 2*P2M * (s) 		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][6] = 2*(r) +1	+ 2*P2M * (s) 		+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][7] = 2*(r) +1	+ 2*P2M * (s)+P2M	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][5] = 2*(r) +1	+ 2*P2M * (s) 		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][8] = 2*(r) +1	+ 2*P2M * (s)+P2M 	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][9] = 2*(r+1)		+ 2*P2M * (s)+P2M 	+ 2*P2M*P2M * (t+1);
                    
                    counter++;
                    
                    (*elementsVec)[counter][0] = 2*(r+1)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][1] = 2*(r)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][2] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][3] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    
                    (*elementsVec)[counter][4] = 2*(r) +1	+ 2*P2M * (s)		+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][6] = 2*(r+1)		+ 2*P2M * (s)+P2M	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][7] = 2*(r+1)		+ 2*P2M * (s)+P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][5] = 2*(r) +1	+ 2*P2M * (s)+P2M	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][8] = 2*(r) +1	+ 2*P2M * (s)+P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][9] = 2*(r+1)		+ 2*P2M * (s+1)		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    
                    counter++;
                    
                    (*elementsVec)[counter][0] = 2*(r)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][1] = 2*(r)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][2] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][3] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    
                    (*elementsVec)[counter][4] = 2*(r) 		+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][6] = 2*(r) +1	+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][7] = 2*(r) +1	+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][5] = 2*(r) +1	+ 2*P2M * (s+1)		+ 2*P2M*P2M * (t);
                    (*elementsVec)[counter][8] = 2*(r) +1	+ 2*P2M * (s+1)		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][9] = 2*(r+1)		+ 2*P2M * (s+1)		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    
                    counter++;
                    
                    (*elementsVec)[counter][0] = 2*(r)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][1] = 2*(r)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][2] = 2*(r)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][3] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    
                    (*elementsVec)[counter][4] = 2*(r)		+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][6] = 2*(r)		+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][7] = 2*(r) +1	+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][5] = 2*(r)		+ 2*P2M*(s+1) 		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][8] = 2*(r) +1	+ 2*P2M*(s+1) 		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][9] = 2*(r) +1 	+ 2*P2M*(s+1) 		+ 2*P2M*P2M * (t+1) ;
                    
                    counter++;
                    
                    (*elementsVec)[counter][0] = 2*(r)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][1] = 2*(r)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][2] = 2*(r)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][3] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    
                    (*elementsVec)[counter][4] = 2*(r)		+ 2*P2M * (s)		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][6] = 2*(r)		+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][7] = 2*(r) +1	+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][5] = 2*(r) 		+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][8] = 2*(r) +1 	+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][9] = 2*(r) +1	+ 2*P2M*(s+1)		+ 2*P2M*P2M * (t+1) ;
                    
                    counter++;
                    
                }
            }
        }
        buildElementsClass(elementsVec, elementFlag);
    }
    else if(FEType == "P1-disc" || FEType == "P1-disc-global")
        buildP1_Disc_Q2_3DCube( N, MM, numProcsCoarseSolve, underlyingLib );
    else if(FEType == "Q1"){
        build3DQ1Cube( N, M, numProcsCoarseSolve, underlyingLib );
    }
    else if(FEType == "Q2"){
        build3DQ2Cube( N, MM, numProcsCoarseSolve, underlyingLib );
    }
    else if(FEType == "Q2-20"){
        build3DQ2_20Cube( N, MM, numProcsCoarseSolve, underlyingLib );
    }
    
    
    
};
    
template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::buildP1_Disc_Q2_3DCube(int N,
                                                        int M,
                                                        int numProcsCoarseSolve,
                                                        std::string underlyingLib){
    
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ScalarTraits;

    
    bool verbose (this->comm_->getRank() == 0);
    
    setRankRange( numProcsCoarseSolve );
    
    if (verbose)
        cout << endl;
    
    SC eps = ScalarTraits<SC>::eps();
    
    int         rank = this->comm_->getRank();
    int         size = this->comm_->getSize() - numProcsCoarseSolve;
    
    SC      h = length/(M*N);//add variable length/width/heigth
    SC      H = length/N;
    
    LO 	nmbPoints_oneDir;
    
    LO nmbElements;
    LO nmbPoints = 4*M*M*M; // 4 points for each element
//        nmbPoints_oneDir 	= N * (M+1) - (N-1);
    
    GO nmbPGlob_oneDir = N * (M+1) - (N-1);
    
    this->numElementsGlob_ = (nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1);
    
    if (rank>=size) {
        M = -1; // keine Schleife wird ausgefuehrt
        nmbElements = 0;
        nmbPoints = 0;
    }
    else{
        nmbElements = M*M*M;
    }
    
    
    int counter = 0;
    int offset_x = (rank % N);
    int offset_y = 0;
    int offset_z = 0;
    
    if ((rank % (N*N))>=N) {
        offset_y = (int) (rank % (N*N))/(N);
    }
    
    if ((rank % (N*N*N))>=N*N ) {
        offset_z = (int) (rank % (N*N*N))/(N*(N));
    }
    
    
    this->pointsRep_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(3,0.0)));
    this->pointsUni_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(3,0.0)));
    this->bcFlagRep_.reset(new std::vector<int> (nmbPoints,0));
    this->bcFlagUni_.reset(new std::vector<int> (nmbPoints,0));
    Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);

    
    if (verbose)
        cout << "-- Building P1-disc Points and Elements according to Q2 ... " << flush;
    
    vec2D_int_ptr_Type elementsVec = Teuchos::rcp(new vec2D_int_Type(nmbElements, vec_int_Type(4, -1)));

    counter = 0;
    LO pCounter = 0;
    GO globalCounterPoints = rank * nmbPoints;
    for (int t=0; t < M; t++) {
        for (int s=0; s < M; s++) {
            for (int r=0; r < M; r++) {

                //point 1
                (*this->pointsRep_)[pCounter][0] = r*h + offset_x * H;
                (*this->pointsRep_)[pCounter][1] = s*h + offset_y * H;
                (*this->pointsRep_)[pCounter][2] = t*h + offset_z * H;
                if ((*this->pointsRep_)[pCounter][0] > (coorRec[0]+length-eps) 	|| (*this->pointsRep_)[pCounter][0] < (coorRec[0]+eps) ||
                    (*this->pointsRep_)[pCounter][1] > (coorRec[1]+width-eps) 	|| (*this->pointsRep_)[pCounter][1] < (coorRec[1]+eps) ||
                    (*this->pointsRep_)[pCounter][2] > (coorRec[2]+height-eps) 	|| (*this->pointsRep_)[pCounter][2] < (coorRec[2]+eps) ) {
                    
                    (*this->bcFlagRep_)[counter] = 1;
                }
                
                (*this->pointsUni_)[pCounter][0] = (*this->pointsRep_)[pCounter][0];
                (*this->pointsUni_)[pCounter][1] = (*this->pointsRep_)[pCounter][1];
                (*this->pointsUni_)[pCounter][2] = (*this->pointsRep_)[pCounter][2];
                
                (*this->bcFlagUni_)[pCounter] = (*this->bcFlagRep_)[pCounter];
                
                (*elementsVec)[counter][0] = pCounter;
                pointsRepGlobMapping[pCounter] = globalCounterPoints;
                globalCounterPoints++;
                pCounter++;
                
                //point 2
                (*this->pointsRep_)[pCounter][0] = (r+1)*h + offset_x * H;
                (*this->pointsRep_)[pCounter][1] = s*h + offset_y * H;
                (*this->pointsRep_)[pCounter][2] = t*h + offset_z * H;
                if ((*this->pointsRep_)[pCounter][0] > (coorRec[0]+length-eps) 	|| (*this->pointsRep_)[pCounter][0] < (coorRec[0]+eps) ||
                    (*this->pointsRep_)[pCounter][1] > (coorRec[1]+width-eps) 	|| (*this->pointsRep_)[pCounter][1] < (coorRec[1]+eps) ||
                    (*this->pointsRep_)[pCounter][2] > (coorRec[2]+height-eps) 	|| (*this->pointsRep_)[pCounter][2] < (coorRec[2]+eps) ) {
                    
                    (*this->bcFlagRep_)[counter] = 1;
                }
                
                (*this->pointsUni_)[pCounter][0] = (*this->pointsRep_)[pCounter][0];
                (*this->pointsUni_)[pCounter][1] = (*this->pointsRep_)[pCounter][1];
                (*this->pointsUni_)[pCounter][2] = (*this->pointsRep_)[pCounter][2];
                
                (*this->bcFlagUni_)[pCounter] = (*this->bcFlagRep_)[pCounter];
                
                (*elementsVec)[counter][1] = pCounter;
                pointsRepGlobMapping[pCounter] = globalCounterPoints;
                globalCounterPoints++;
                pCounter++;
                //point 3
                (*this->pointsRep_)[pCounter][0] = r*h + offset_x * H;
                (*this->pointsRep_)[pCounter][1] = (s+1)*h + offset_y * H;
                (*this->pointsRep_)[pCounter][2] = t*h + offset_z * H;
                if ((*this->pointsRep_)[pCounter][0] > (coorRec[0]+length-eps) 	|| (*this->pointsRep_)[pCounter][0] < (coorRec[0]+eps) ||
                    (*this->pointsRep_)[pCounter][1] > (coorRec[1]+width-eps) 	|| (*this->pointsRep_)[pCounter][1] < (coorRec[1]+eps) ||
                    (*this->pointsRep_)[pCounter][2] > (coorRec[2]+height-eps) 	|| (*this->pointsRep_)[pCounter][2] < (coorRec[2]+eps) ) {
                    
                    (*this->bcFlagRep_)[counter] = 1;
                }
                
                (*this->pointsUni_)[pCounter][0] = (*this->pointsRep_)[pCounter][0];
                (*this->pointsUni_)[pCounter][1] = (*this->pointsRep_)[pCounter][1];
                (*this->pointsUni_)[pCounter][2] = (*this->pointsRep_)[pCounter][2];
                
                (*this->bcFlagUni_)[pCounter] = (*this->bcFlagRep_)[pCounter];

                (*elementsVec)[counter][2] = pCounter;
                pointsRepGlobMapping[pCounter] = globalCounterPoints;
                globalCounterPoints++;
                pCounter++;

                //point 4
                (*this->pointsRep_)[pCounter][0] = r*h + offset_x * H;
                (*this->pointsRep_)[pCounter][1] = s*h + offset_y * H;
                (*this->pointsRep_)[pCounter][2] = (t+1)*h + offset_z * H;
                if ((*this->pointsRep_)[pCounter][0] > (coorRec[0]+length-eps) 	|| (*this->pointsRep_)[pCounter][0] < (coorRec[0]+eps) ||
                    (*this->pointsRep_)[pCounter][1] > (coorRec[1]+width-eps) 	|| (*this->pointsRep_)[pCounter][1] < (coorRec[1]+eps) ||
                    (*this->pointsRep_)[pCounter][2] > (coorRec[2]+height-eps) 	|| (*this->pointsRep_)[pCounter][2] < (coorRec[2]+eps) ) {
                    
                    (*this->bcFlagRep_)[counter] = 1;
                }
                
                (*this->pointsUni_)[pCounter][0] = (*this->pointsRep_)[pCounter][0];
                (*this->pointsUni_)[pCounter][1] = (*this->pointsRep_)[pCounter][1];
                (*this->pointsUni_)[pCounter][2] = (*this->pointsRep_)[pCounter][2];
                
                (*this->bcFlagUni_)[pCounter] = (*this->bcFlagRep_)[pCounter];
                
                (*elementsVec)[counter][3] = pCounter;
                pointsRepGlobMapping[pCounter] = globalCounterPoints;
                globalCounterPoints++;
                pCounter++;
                
                counter++;
            }
        }
    }
    this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );
    
    this->mapUnique_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );

    buildElementsClass(elementsVec);
    
    if (verbose)
        cout << "done!" << endl;

}
    
template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::build3DQ1Cube(int N,
                                                int M,
                                                int numProcsCoarseSolve,
                                                std::string underlyingLib)
{
    
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ScalarTraits;
    
    bool verbose (this->comm_->getRank() == 0);
    
    if (verbose)
        cout << endl;
    
    setRankRange( numProcsCoarseSolve );
    
    SC eps = ScalarTraits<SC>::eps();
    
    int         rank = this->comm_->getRank();
    int         size = this->comm_->getSize() - numProcsCoarseSolve;
    
    SC      h = length/(M*N);//add variable length/width/heigth
    SC      H = length/N;
    
    LO nmbElements;

    LO nmbPoints_oneDir 	= N * (M+1) - (N-1);
    LO nmbPoints			= (M+1)*(M+1)*(M+1);

    GO nmbPGlob_oneDir = N * (M+1) - (N-1);
    
    this->numElementsGlob_ = (nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1);
    
    if (rank>=size) {
        M = -1; // keine Schleife wird ausgefuehrt
        nmbElements = 0;
        nmbPoints = 0;
    }
    else{
        nmbElements = M*M*M;
    }
    
    this->pointsRep_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(3,0.0)));
    this->bcFlagRep_.reset(new std::vector<int> (nmbPoints,0));
    vec2D_int_ptr_Type elementsVec = Teuchos::rcp(new std::vector<std::vector<int> >(nmbElements, std::vector<int>(8, -1)));
    Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);
    
    int counter = 0;
    int offset_x = (rank % N);
    int offset_y = 0;
    int offset_z = 0;
    
    if ((rank % (N*N))>=N) {
        offset_y = (int) (rank % (N*N))/(N);
    }
    
    if ((rank % (N*N*N))>=N*N ) {
        offset_z = (int) (rank % (N*N*N))/(N*(N));
    }
    for (int t=0; t < M+1; t++) {
        for (int s=0; s < M+1; s++) {
            for (int r=0; r < M+1; r++) {
                (*this->pointsRep_)[counter][0] = r*h + offset_x * H;
                (*this->pointsRep_)[counter][1] = s*h + offset_y * H;
                (*this->pointsRep_)[counter][2] = t*h + offset_z * H;
                
                pointsRepGlobMapping[counter] = r + s*nmbPoints_oneDir + t*nmbPoints_oneDir*nmbPoints_oneDir \
                + offset_x*(M) + offset_y*(nmbPoints_oneDir)*M + offset_z*(nmbPoints_oneDir)*(nmbPoints_oneDir)*M  ;
                
                if ((*this->pointsRep_)[counter][0] > (coorRec[0]+length-eps) 	|| (*this->pointsRep_)[counter][0] < (coorRec[0]+eps) ||
                    (*this->pointsRep_)[counter][1] > (coorRec[1]+width-eps) 	|| (*this->pointsRep_)[counter][1] < (coorRec[1]+eps) ||
                    (*this->pointsRep_)[counter][2] > (coorRec[2]+height-eps) 	|| (*this->pointsRep_)[counter][2] < (coorRec[2]+eps) ) {
                    
                    (*this->bcFlagRep_)[counter] = 1;
                }
                counter++;
            }
        }
    }
    this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );
    
    this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );
    
    this->pointsUni_.reset(new std::vector<std::vector<double> >(this->mapUnique_->getNodeNumElements(),std::vector<double>(3,0.0)));
    this->bcFlagUni_.reset(new std::vector<int> (this->mapUnique_->getNodeNumElements(),0));
    
    for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
        
        LO index = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement(i) );
        
        (*this->pointsUni_)[i][0] = (*this->pointsRep_)[index][0];
        (*this->pointsUni_)[i][1] = (*this->pointsRep_)[index][1];
        (*this->pointsUni_)[i][2] = (*this->pointsRep_)[index][2];
        (*this->bcFlagUni_)[i] = (*this->bcFlagRep_)[index];
    }
    
    LO offset = (M+1);
    
    counter = 0;
    for (int t=0; t < M; t++) {
        for (int s=0; s < M; s++) {
            for (int r=0; r < M; r++) {
                
                (*elementsVec)[counter][0] = r      + (M+1) * (s)	+ (M+1)*(M+1) * t ;
                (*elementsVec)[counter][1] = r + 1  + (M+1) * (s)	+ (M+1)*(M+1) * t ;
                (*elementsVec)[counter][2] = r + 1  + (M+1) * (s+1)	+ (M+1)*(M+1) * t ;
                (*elementsVec)[counter][3] = r      + (M+1) * (s+1)	+ (M+1)*(M+1) * t ;
                
                (*elementsVec)[counter][4] = r      + (M+1) * (s)	+ (M+1)*(M+1) * (t+1) ;
                (*elementsVec)[counter][5] = r + 1  + (M+1) * (s)	+ (M+1)*(M+1) * (t+1) ;
                (*elementsVec)[counter][6] = r + 1  + (M+1) * (s+1)	+ (M+1)*(M+1) * (t+1) ;
                (*elementsVec)[counter][7] = r      + (M+1) * (s+1)	+ (M+1)*(M+1) * (t+1) ;

                counter++;
                
            }
        }
    }
    buildElementsClass(elementsVec);
}

    
template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::build3DQ2Cube(int N,
                                                int M,
                                                int numProcsCoarseSolve,
                                                std::string underlyingLib)
{
   
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ScalarTraits;
    
    bool verbose (this->comm_->getRank() == 0);
    
    if (verbose)
        cout << endl;
    
    setRankRange( numProcsCoarseSolve );
    
    SC eps = ScalarTraits<SC>::eps();
    
    int         rank = this->comm_->getRank();
    int         size = this->comm_->getSize() - numProcsCoarseSolve;
    
    SC      h = length/(M*N);//add variable length/width/heigth
    SC      H = length/N;
    
    LO nmbElements;

    LO nmbPoints_oneDir =  N * (2*(M+1)-1) - (N-1);
    LO nmbPoints = (2*(M+1)-1)*(2*(M+1)-1)*(2*(M+1)-1);
    
    GO nmbPGlob_oneDir = N * (M+1) - (N-1);
    
    this->numElementsGlob_ = (nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1);
    
    if (rank>=size) {
        M = -1; // keine Schleife wird ausgefuehrt
        nmbElements = 0;
        nmbPoints = 0;
    }
    else{
        nmbElements = M*M*M;
    }
    
    this->pointsRep_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(3,0.0)));
    this->bcFlagRep_.reset(new std::vector<int> (nmbPoints,0));
    vec2D_int_ptr_Type elementsVec = Teuchos::rcp(new std::vector<std::vector<int> >(nmbElements,std::vector<int>(27,-1)));
    Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);
    
    int counter = 0;
    int offset_x = (rank % N);
    int offset_y = 0;
    int offset_z = 0;
    
    if ((rank % (N*N))>=N) {
        offset_y = (int) (rank % (N*N))/(N);
    }
    
    if ((rank % (N*N*N))>=N*N ) {
        offset_z = (int) (rank % (N*N*N))/(N*(N));
    }
    bool p1point;
    int p1_s = 0;
    int p1_r = 0;
    int p1_t = 0;
    for (int t=0; t < 2*(M+1)-1; t++) {
        for (int s=0; s < 2*(M+1)-1; s++) {
            for (int r=0; r < 2*(M+1)-1; r++) {
                p1point = false;
                if (s%2==0 && r%2==0 && t%2==0) {
                    p1point = true;
                    p1_s = s/2;
                    p1_r = r/2;
                    p1_t = t/2;
                }
                (*this->pointsRep_)[counter][0] = r*h/2.0 + offset_x * H;
                if ((*this->pointsRep_)[counter][0]<eps && (*this->pointsRep_)[counter][0]>-eps) (*this->pointsRep_)[counter][0]=0.0;
                (*this->pointsRep_)[counter][1] = s*h/2.0 + offset_y * H;
                if ((*this->pointsRep_)[counter][1]<eps && (*this->pointsRep_)[counter][1]>-eps) (*this->pointsRep_)[counter][1]=0.0;
                (*this->pointsRep_)[counter][2] = t*h/2.0 + offset_z * H;
                if ((*this->pointsRep_)[counter][2]<eps && (*this->pointsRep_)[counter][2]>-eps) (*this->pointsRep_)[counter][2]=0.0;
                
                pointsRepGlobMapping[counter] = r + s*nmbPoints_oneDir + t*nmbPoints_oneDir*nmbPoints_oneDir \
                + offset_x*(2*(M+1)-2) + offset_y*(nmbPoints_oneDir)*(2*(M+1)-2) + offset_z*(nmbPoints_oneDir)*(nmbPoints_oneDir)*(2*(M+1)-2) ;
                
                if ((*this->pointsRep_)[counter][0] > (coorRec[0]+length-eps) || (*this->pointsRep_)[counter][0] < (coorRec[0]+eps) ||
                    (*this->pointsRep_)[counter][1] > (coorRec[1]+width-eps) 	|| (*this->pointsRep_)[counter][1] < (coorRec[1]+eps) ||
                    (*this->pointsRep_)[counter][2] > (coorRec[2]+height-eps) || (*this->pointsRep_)[counter][2] < (coorRec[2]+eps) ) {
                    (*this->bcFlagRep_)[counter] = 1;
                    
                }
                
                counter++;
            }
        }
    }
    
    this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );
    
    this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );
    
    this->pointsUni_.reset(new std::vector<std::vector<double> >(this->mapUnique_->getNodeNumElements(),std::vector<double>(3,0.0)));
    this->bcFlagUni_.reset(new std::vector<int> (this->mapUnique_->getNodeNumElements(),0));
    
    for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
        (*this->pointsUni_)[i][0] = (this->mapUnique_->getGlobalElement(i) % nmbPoints_oneDir) * h/2;
        if ((*this->pointsUni_)[i][0]<eps && (*this->pointsUni_)[i][0]>-eps) (*this->pointsUni_)[i][0]=0.0;
        
        (*this->pointsUni_)[i][1] = ((int) ((this->mapUnique_->getGlobalElement(i) % (nmbPoints_oneDir*nmbPoints_oneDir)) / nmbPoints_oneDir) + eps) *h/2;
        if ((*this->pointsUni_)[i][1]<eps && (*this->pointsUni_)[i][1]>-eps) (*this->pointsUni_)[i][1]=0.0;
        
        (*this->pointsUni_)[i][2] = ((int)(this->mapUnique_->getGlobalElement(i) / (nmbPoints_oneDir*nmbPoints_oneDir) + eps)) * h/2;
        if ((*this->pointsUni_)[i][2]<eps && (*this->pointsUni_)[i][2]>-eps) (*this->pointsUni_)[i][2]=0.0;
        
        if ((*this->pointsUni_)[i][0] > (coorRec[0]+length-eps) 	|| (*this->pointsUni_)[i][0] < (coorRec[0]+eps) ||
            (*this->pointsUni_)[i][1] > (coorRec[1]+width-eps) 	|| (*this->pointsUni_)[i][1] < (coorRec[1]+eps) ||
            (*this->pointsUni_)[i][2] > (coorRec[2]+height-eps) 	|| (*this->pointsUni_)[i][2] < (coorRec[2]+eps) ) {
            (*this->bcFlagUni_)[i] = 1;
            
        }
    }
    
    int    P2M = 2*(M+1)-1;
    
    counter = 0;
    for (int t=0; t < M; t++) {
        for (int s=0; s < M; s++) {
            for (int r=0; r < M; r++) {
                
                (*elementsVec)[counter][0] = 2*(r)      + 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                (*elementsVec)[counter][1] = 2*(r+1)    + 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                (*elementsVec)[counter][2] = 2*(r+1)    + 2*P2M * (s+1)	+ 2*P2M*P2M * (t) ;
                (*elementsVec)[counter][3] = 2*(r)      + 2*P2M * (s+1)	+ 2*P2M*P2M * (t) ;
                
                (*elementsVec)[counter][4] = 2*(r)      + 2*P2M * (s)	+ 2*P2M*P2M * (t+1) ;
                (*elementsVec)[counter][5] = 2*(r+1)    + 2*P2M * (s)	+ 2*P2M*P2M * (t+1) ;
                (*elementsVec)[counter][6] = 2*(r+1)    + 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                (*elementsVec)[counter][7] = 2*(r)      + 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                
                (*elementsVec)[counter][8] = 2*(r)+1    + 2*P2M * (s)           + 2*P2M*P2M * (t) ;
                (*elementsVec)[counter][9] = 2*(r+1)    + 2*P2M * (s)   + P2M	+ 2*P2M*P2M * (t) ;
                (*elementsVec)[counter][10] = 2*(r)+1    + 2*P2M * (s+1)        + 2*P2M*P2M * (t) ;
                (*elementsVec)[counter][11] = 2*(r)      + 2*P2M * (s)  + P2M	+ 2*P2M*P2M * (t) ;

                (*elementsVec)[counter][12] = 2*(r)+1    + 2*P2M * (s)           + 2*P2M*P2M * (t+1) ;
                (*elementsVec)[counter][13] = 2*(r+1)    + 2*P2M * (s)   + P2M	+ 2*P2M*P2M * (t+1) ;
                (*elementsVec)[counter][14] = 2*(r)+1    + 2*P2M * (s+1)        + 2*P2M*P2M * (t+1) ;
                (*elementsVec)[counter][15] = 2*(r)      + 2*P2M * (s)  + P2M	+ 2*P2M*P2M * (t+1) ;

                (*elementsVec)[counter][16] = 2*(r)      + 2*P2M * (s)	+ 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][17] = 2*(r+1)    + 2*P2M * (s)	+ 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][18] = 2*(r+1)    + 2*P2M * (s+1)	+ 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][19] = 2*(r)      + 2*P2M * (s+1)	+ 2*P2M*P2M * (t) + P2M*P2M;

                (*elementsVec)[counter][20] = 2*(r)+1    + 2*P2M * (s)           + 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][21] = 2*(r+1)    + 2*P2M * (s)   + P2M	+ 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][22] = 2*(r)+1    + 2*P2M * (s+1)        + 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][23] = 2*(r)      + 2*P2M * (s)  + P2M	+ 2*P2M*P2M * (t) + P2M*P2M;
                
                (*elementsVec)[counter][24] = 2*(r)+1    + 2*P2M * (s)   + P2M  + 2*P2M*P2M * (t);
                (*elementsVec)[counter][25] = 2*(r)+1    + 2*P2M * (s)   + P2M  + 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][26] = 2*(r)+1    + 2*P2M * (s)   + P2M  + 2*P2M*P2M * (t+1);

                counter++;
            }
        }
    }
    buildElementsClass(elementsVec);
}
    
template <class SC, class LO, class GO, class NO>
GO MeshStructured<SC,LO,GO,NO>::globalID_Q2_20Cube(int r, int s , int t, int &rr, int off_x, int off_y, int off_z, int M, int N, GO nmbPoints_oneDirFull, GO nmbPoints_oneDirMid){

    GO index = -1;
    bool setPoint = false;

    if (r%2==0 && t%2==0)
        setPoint = true;
    else{
        if (s%2==0 && t%2==0)
            setPoint = true;
        else{
            if (t%2==1 && r%2==0 && s%2==0) {
                setPoint = true;
            }
        }
    }
    
    
    if (setPoint) {
        
        long long sizeFullSquare = nmbPoints_oneDirMid * nmbPoints_oneDirFull + nmbPoints_oneDirMid * (nmbPoints_oneDirMid-1);
        long long sizeNotFullSquare = nmbPoints_oneDirMid * nmbPoints_oneDirMid ;
        int ss = s/2;
        int tt = t/2;

        index = rr + nmbPoints_oneDirFull * ( ss + (s%2) ) * (t%2==0) + nmbPoints_oneDirMid * ss;
        index += sizeFullSquare * ( tt +  (t%2) ) + sizeNotFullSquare * tt;
        if (s%2==0)
            index += off_x * 2*M ;
        else
            index += off_x * M ;

        index += off_y * ( nmbPoints_oneDirFull * M + nmbPoints_oneDirMid * M );
        index += off_z * ( sizeFullSquare * M + sizeNotFullSquare * M );
        rr++;
    }
    
    return index;
}

template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::build3DQ2_20Cube(int N,
                                                   int M,
                                                   int numProcsCoarseSolve,
                                                   std::string underlyingLib)
{
    
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ScalarTraits;
    
    bool verbose (this->comm_->getRank() == 0);
    
    setRankRange( numProcsCoarseSolve );
    
    if (verbose)
        cout << endl;
    if (verbose)
        std::cout << "WARNING! Not working properly in parallel - fix global indexing." << std::endl;
    
    SC eps = ScalarTraits<SC>::eps();
    
    int         rank = this->comm_->getRank();
    int         size = this->comm_->getSize() - numProcsCoarseSolve;
    
    SC      h = length/(M*N);//add variable length/width/heigth
    SC      H = length/N;
    
    LO nmbElements;
    
    LO nmbPoints_oneDirFull =  N * (2*(M+1)-1) - (N-1);
    LO nmbPoints_oneDirMid =  N * (M+1) - (N-1);
    
    LO nmbPoints = ( (M+1) * (2*(M+1)-1) + M * (M+1) ) * (M+1) +
                   ( (M+1) * (M+1) ) * M;
    
    GO nmbPGlob_oneDir = N * (M+1) - (N-1);
    
    this->numElementsGlob_ = (nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1);
    
    if (rank>=size) {
        M = -1; // keine Schleife wird ausgefuehrt
        nmbElements = 0;
        nmbPoints = 0;
    }
    else{
        nmbElements = M*M*M;
    }
    
    this->pointsRep_.reset(new std::vector<std::vector<double> >(0));
    this->bcFlagRep_.reset(new std::vector<int> (0));
    vec2D_int_ptr_Type elementsVec = Teuchos::rcp(new std::vector<std::vector<int> >(nmbElements, std::vector<int>(20,-1)));
    Teuchos::Array<GO> pointsRepGlobMapping(0);
    
    int counter = 0;
    int offset_x = (rank % N);
    int offset_y = 0;
    int offset_z = 0;
    
    if ((rank % (N*N))>=N) {
        offset_y = (int) (rank % (N*N))/(N);
    }
    
    if ((rank % (N*N*N))>=N*N ) {
        offset_z = (int) (rank % (N*N*N))/(N*(N));
    }
    int rr=0;
    for (int t=0; t < 2*(M+1)-1; t++) {
        for (int s=0; s < 2*(M+1)-1; s++) {
            rr=0;
            for (int r=0; r < 2*(M+1)-1; r++) {
                GO index = globalID_Q2_20Cube( r, s, t, rr, offset_x, offset_y, offset_z, M, N,
                                              nmbPoints_oneDirFull, nmbPoints_oneDirMid );
                
                if ( index>-1 ) {
                    std::vector<double> p(3,0.0);
                    p[0] = r*h/2.0 + offset_x * H;
                    p[1] = s*h/2.0 + offset_y * H;
                    p[2] = t*h/2.0 + offset_z * H;
                    this->pointsRep_->push_back(p);
                    pointsRepGlobMapping.push_back( index );

                    if (p[0] > (coorRec[0]+length-eps) || p[0] < (coorRec[0]+eps) ||
                        p[1] > (coorRec[1]+width-eps)  || p[1] < (coorRec[1]+eps) ||
                        p[2] > (coorRec[2]+height-eps) || p[2] < (coorRec[2]+eps) )
                        this->bcFlagRep_->push_back(1);
                    else
                        this->bcFlagRep_->push_back(0);
                    counter++;
                }
            }
        }
    }

    this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );

    this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );
    
    this->pointsUni_.reset(new std::vector<std::vector<double> >(this->mapUnique_->getNodeNumElements(),std::vector<double>(3,0.0)));
    this->bcFlagUni_.reset(new std::vector<int> (this->mapUnique_->getNodeNumElements(),0));
    
    for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
        
        LO index = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement(i) );
        
        (*this->pointsUni_)[i][0] = (*this->pointsRep_)[index][0];
        (*this->pointsUni_)[i][1] = (*this->pointsRep_)[index][1];
        (*this->pointsUni_)[i][2] = (*this->pointsRep_)[index][2];
        (*this->bcFlagUni_)[i] = (*this->bcFlagRep_)[index];
    }
    
    int    P2M = 2*(M+1)-1;
    int    P1M = M+1;
    counter = 0;
    for (int t=0; t < M; t++) {
        for (int s=0; s < M; s++) {
            for (int r=0; r < M; r++) {

                (*elementsVec)[counter][0] = 2*(r)      + (P2M+P1M) * (s) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t) ;
                (*elementsVec)[counter][1] = 2*(r+1)    + (P2M+P1M) * (s) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t) ;
                (*elementsVec)[counter][2] = 2*(r+1)    + (P2M+P1M) * (s+1) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t) ;
                (*elementsVec)[counter][3] = 2*(r)      + (P2M+P1M) * (s+1) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t) ;
                
                (*elementsVec)[counter][4] = 2*(r)      + (P2M+P1M) * (s) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t+1) ;
                (*elementsVec)[counter][5] = 2*(r+1)    + (P2M+P1M) * (s) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t+1) ;
                (*elementsVec)[counter][6] = 2*(r+1)    + (P2M+P1M) * (s+1) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t+1) ;
                (*elementsVec)[counter][7] = 2*(r)      + (P2M+P1M) * (s+1) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t+1) ;
                
                (*elementsVec)[counter][8] = 2*(r)+1        + (P2M+P1M) * (s) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t) ;
                (*elementsVec)[counter][9] = (r+1)          + (P2M+P1M) * (s) + P2M + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t) ;
                (*elementsVec)[counter][10] = 2*(r)+1       + (P2M+P1M) * (s+1) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t) ;
                (*elementsVec)[counter][11] = r             + (P2M+P1M) * (s) + P2M + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t) ;
                
                (*elementsVec)[counter][12] = 2*(r)+1        + (P2M+P1M) * (s) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t+1) ;
                (*elementsVec)[counter][13] = (r+1)          + (P2M+P1M) * (s) + P2M + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t+1) ;
                (*elementsVec)[counter][14] = 2*(r)+1       + (P2M+P1M) * (s+1) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t+1) ;
                (*elementsVec)[counter][15] = r             + (P2M+P1M) * (s) + P2M + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t+1) ;
                
                (*elementsVec)[counter][16] = (r)        + (P1M) * (s) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t) + P2M*(M+1)+P1M*M;
                (*elementsVec)[counter][17] = (r+1)      + (P1M) * (s) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t) + P2M*(M+1)+P1M*M;
                (*elementsVec)[counter][18] = (r+1)      + (P1M) * (s+1) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t) + P2M*(M+1)+P1M*M;
                (*elementsVec)[counter][19] = r          + (P1M) * (s+1) + ( (P2M*(M+1)+P1M*M) + (P1M*P1M) ) * (t) + P2M*(M+1)+P1M*M;
                                
                counter++;
            }
        }
    }

}
    
template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::build3DQ2BFS(int N,
                                                int M,
                                                int numProcsCoarseSolve,
                                                std::string underlyingLib){

    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ScalarTraits;
    
    typedef ScalarTraits<SC> ST;
    SC eps = ST::eps();
    
    bool verbose (this->comm_->getRank() == 0);
    
    setRankRange( numProcsCoarseSolve );
    
    int         rank = this->comm_->getRank();
    int         size = this->comm_->getSize() - numProcsCoarseSolve;
    
    int         bfs_multiplier = (int) 2*(length)-1;
    
    int         nmbSubdomainsSquares = size / bfs_multiplier;
    int         nmbSubdomainsSquares_OneDir = (std::pow(nmbSubdomainsSquares,1./3.) + 100*eps); // same as N
    
    SC      h = ST::one()/(M*N);
    SC      H = ST::one()/N;
    

    LO nmbPoints_oneDir 	=  N * (2*(M+1)-1) - (N-1);
    LO nmbPoints_oneDir_allSubdomain = length * nmbPoints_oneDir - (length-1) ;
    LO  nmbPoints 			= (2*(M+1)-1)*(2*(M+1)-1)*(2*(M+1)-1);
    
    
    GO nmbPGlob_oneDir = N * (M+1) - (N-1);
    this->numElementsGlob_ = (nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1) * bfs_multiplier;
    LO nmbElements;
    if (rank>=size) {
        M = -1; // keine Schleife wird ausgefuehrt
        nmbElements = 0;
        nmbPoints = 0;
    }
    else{
        nmbElements = M*M*M;
    }
    
    this->pointsRep_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(3,0.0)));
    this->bcFlagRep_.reset(new std::vector<int> (nmbPoints,0));
    vec2D_int_ptr_Type elementsVec = Teuchos::rcp(new std::vector<std::vector<int> >(nmbElements, std::vector<int>(27,-1)));
    Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);
    
    int whichSquareSet = (int)rank / nmbSubdomainsSquares;
    
    int offset_Squares_x = (int) (whichSquareSet+1) / 2;
    int offset_Squares_y = 0;
    int offset_Squares_z = ((whichSquareSet+1) % 2);
    
    int counter = 0;
    int offset_x = ((rank - nmbSubdomainsSquares*whichSquareSet) % N);
    int offset_y = 0;
    int offset_z = 0;
    if (((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))>=N) {
        offset_y = (int) ((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))/(N);
    }
    
    if (((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N*N))>=N*N ) {
        offset_z = (int) ((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N*N))/(N*(N));
    }
    
    for (int t=0; t < 2*(M+1)-1; t++) {
        for (int s=0; s < 2*(M+1)-1; s++) {
            for (int r=0; r < 2*(M+1)-1; r++) {
                
                (*this->pointsRep_)[counter][0] = coorRec[0] + r*h/2. + offset_x * H + offset_Squares_x * H * nmbSubdomainsSquares_OneDir;
                
                (*this->pointsRep_)[counter][1] = coorRec[1] + s*h/2. + offset_y * H + offset_Squares_y * H * nmbSubdomainsSquares_OneDir;
                
                (*this->pointsRep_)[counter][2] = coorRec[2] + t*h/2. + offset_z * H + offset_Squares_z * H * nmbSubdomainsSquares_OneDir;
                
                pointsRepGlobMapping[counter] = r
                + s*(nmbPoints_oneDir_allSubdomain);
                if (offset_Squares_x > 0 && offset_Squares_z == 0 && offset_z+1!=nmbSubdomainsSquares_OneDir) {
                    pointsRepGlobMapping[counter] -= s*(nmbPoints_oneDir-1);
                }
                else if(offset_Squares_x > 0 && offset_Squares_z == 0 && offset_z+1==nmbSubdomainsSquares_OneDir && t!=2*M){
                    pointsRepGlobMapping[counter] -= s*(nmbPoints_oneDir-1);
                }
                
                pointsRepGlobMapping[counter] += t*nmbPoints_oneDir_allSubdomain*nmbPoints_oneDir;
                if (offset_Squares_x > 0 && offset_Squares_z == 0 ) {
                    pointsRepGlobMapping[counter] -= t*(nmbPoints_oneDir-1)*(nmbPoints_oneDir);
                }
                
                pointsRepGlobMapping[counter] += offset_x*(2*M)
                + offset_y*( nmbPoints_oneDir_allSubdomain * 2*M );
                if (offset_Squares_x > 0 && offset_Squares_z == 0 && offset_z+1!=nmbSubdomainsSquares_OneDir) {
                    pointsRepGlobMapping[counter] -= offset_y*2*M*(nmbPoints_oneDir-1);
                }
                else if(offset_Squares_x > 0 && offset_Squares_z == 0 && offset_z+1==nmbSubdomainsSquares_OneDir && t!=2*M){
                    pointsRepGlobMapping[counter] -= offset_y*2*M*(nmbPoints_oneDir-1);
                }
                
                pointsRepGlobMapping[counter] += offset_z * 2*M * nmbPoints_oneDir_allSubdomain * nmbPoints_oneDir;
                if (offset_Squares_x > 0 && offset_Squares_z == 0 ) {
                    pointsRepGlobMapping[counter] -= offset_z*2*M*(nmbPoints_oneDir-1)*(nmbPoints_oneDir);
                }
                
                pointsRepGlobMapping[counter] += offset_Squares_x * 2*M * nmbSubdomainsSquares_OneDir;
                if (offset_Squares_z == 0 && offset_z+1!=nmbSubdomainsSquares_OneDir) {
                    pointsRepGlobMapping[counter] -= 2*M * nmbSubdomainsSquares_OneDir;
                }
                else if(offset_Squares_z == 0 && offset_z+1==nmbSubdomainsSquares_OneDir && t!=2*M){
                    pointsRepGlobMapping[counter] -= 2*M * nmbSubdomainsSquares_OneDir;
                }
                
                pointsRepGlobMapping[counter] += offset_Squares_z * nmbPoints_oneDir_allSubdomain * ((2*M) * nmbSubdomainsSquares_OneDir+1) * 2*M * nmbSubdomainsSquares_OneDir;
                if (offset_Squares_z > 0 ) {
                    pointsRepGlobMapping[counter] -= (nmbPoints_oneDir-1)*(nmbPoints_oneDir-1)*(nmbPoints_oneDir);
                }
                counter++;
            }
        }
    }
    
    this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );
    
    this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );
    
    if (verbose)
        cout << "-- Building Q2 Unique Points ... " << flush;
    
    this->pointsUni_.reset(new std::vector<std::vector<double> >(this->mapUnique_->getNodeNumElements(),std::vector<double>(3,0.0)));
    this->bcFlagUni_.reset(new std::vector<int> (this->mapUnique_->getNodeNumElements(),0));
    
    LO index;
    for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
        
        index = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement(i) );
        
        (*this->pointsUni_)[i][0] = (*this->pointsRep_)[index][0];
        (*this->pointsUni_)[i][1] = (*this->pointsRep_)[index][1];
        (*this->pointsUni_)[i][2] = (*this->pointsRep_)[index][2];
        (*this->bcFlagUni_)[i] = (*this->bcFlagRep_)[index];
    }
    
    if (verbose)
        cout << " done! --" << endl;
    
    int    P2M = 2*(M+1)-1;
    
    counter = 0;
    for (int t=0; t < M; t++) {
        for (int s=0; s < M; s++) {
            for (int r=0; r < M; r++) {
                
                (*elementsVec)[counter][0] = 2*(r)      + 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                (*elementsVec)[counter][1] = 2*(r+1)    + 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                (*elementsVec)[counter][2] = 2*(r+1)    + 2*P2M * (s+1)	+ 2*P2M*P2M * (t) ;
                (*elementsVec)[counter][3] = 2*(r)      + 2*P2M * (s+1)	+ 2*P2M*P2M * (t) ;
                
                (*elementsVec)[counter][4] = 2*(r)      + 2*P2M * (s)	+ 2*P2M*P2M * (t+1) ;
                (*elementsVec)[counter][5] = 2*(r+1)    + 2*P2M * (s)	+ 2*P2M*P2M * (t+1) ;
                (*elementsVec)[counter][6] = 2*(r+1)    + 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                (*elementsVec)[counter][7] = 2*(r)      + 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                
                (*elementsVec)[counter][8] = 2*(r)+1    + 2*P2M * (s)           + 2*P2M*P2M * (t) ;
                (*elementsVec)[counter][9] = 2*(r+1)    + 2*P2M * (s)   + P2M	+ 2*P2M*P2M * (t) ;
                (*elementsVec)[counter][10] = 2*(r)+1    + 2*P2M * (s+1)        + 2*P2M*P2M * (t) ;
                (*elementsVec)[counter][11] = 2*(r)      + 2*P2M * (s)  + P2M	+ 2*P2M*P2M * (t) ;
                
                (*elementsVec)[counter][12] = 2*(r)+1    + 2*P2M * (s)           + 2*P2M*P2M * (t+1) ;
                (*elementsVec)[counter][13] = 2*(r+1)    + 2*P2M * (s)   + P2M	+ 2*P2M*P2M * (t+1) ;
                (*elementsVec)[counter][14] = 2*(r)+1    + 2*P2M * (s+1)        + 2*P2M*P2M * (t+1) ;
                (*elementsVec)[counter][15] = 2*(r)      + 2*P2M * (s)  + P2M	+ 2*P2M*P2M * (t+1) ;
                
                (*elementsVec)[counter][16] = 2*(r)      + 2*P2M * (s)	+ 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][17] = 2*(r+1)    + 2*P2M * (s)	+ 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][18] = 2*(r+1)    + 2*P2M * (s+1)	+ 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][19] = 2*(r)      + 2*P2M * (s+1)	+ 2*P2M*P2M * (t) + P2M*P2M;
                
                (*elementsVec)[counter][20] = 2*(r)+1    + 2*P2M * (s)           + 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][21] = 2*(r+1)    + 2*P2M * (s)   + P2M	+ 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][22] = 2*(r)+1    + 2*P2M * (s+1)        + 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][23] = 2*(r)      + 2*P2M * (s)  + P2M	+ 2*P2M*P2M * (t) + P2M*P2M;
                
                (*elementsVec)[counter][24] = 2*(r)+1    + 2*P2M * (s)   + P2M  + 2*P2M*P2M * (t);
                (*elementsVec)[counter][25] = 2*(r)+1    + 2*P2M * (s)   + P2M  + 2*P2M*P2M * (t) + P2M*P2M;
                (*elementsVec)[counter][26] = 2*(r)+1    + 2*P2M * (s)   + P2M  + 2*P2M*P2M * (t+1);
                
                counter++;
                
            }
        }
    }
}

template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::buildMesh2DBFS(std::string FEType,
                                                    int N,
                                                    int M,
                                                    int numProcsCoarseSolve,
                                                    std::string underlyingLib) {

    
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ScalarTraits;
    
    TEUCHOS_TEST_FOR_EXCEPTION(!(M>=1),std::logic_error,"H/h is to small.");
    TEUCHOS_TEST_FOR_EXCEPTION(this->comm_.is_null(),std::runtime_error,"comm_ is null.");
    
    SC eps = ScalarTraits<SC>::eps();
    
    bool verbose (this->comm_->getRank() == 0);
    
    setRankRange( numProcsCoarseSolve );
    
    int         rank = this->comm_->getRank();
    int         size = this->comm_->getSize() - numProcsCoarseSolve;
    
    int         bfs_multiplier = (int) 2*(length)-1;

    int         nmbSubdomainsSquares = size / bfs_multiplier;
    int         nmbSubdomainsSquares_OneDir = (std::pow(nmbSubdomainsSquares,1/2.) + 100*eps);
    
    vec2D_int_ptr_Type elementsVec;
    
    LO nmbElements;
    LO nmbPoints;
    
    double      h = 1./(M*N);
    double      H = 1./N;
    GO 	nmbPoints_oneDir;
    GO  nmbPoints_oneDir_allSubdomain;
    if (FEType == "P0") {
        nmbPoints_oneDir 	= N * (M+1) - (N-1) ;
        nmbPoints_oneDir_allSubdomain 	= length * nmbPoints_oneDir - (length-1) ;
        nmbPoints			= (M+1)*(M+1) ;
    }
    else if (FEType == "P1") {
        nmbPoints_oneDir 	= N * (M+1) - (N-1) ;
        nmbPoints_oneDir_allSubdomain 	= length * nmbPoints_oneDir - (length-1) ;
        nmbPoints			= (M+1)*(M+1) ;
    }
    else if(FEType == "P2"){
        nmbPoints_oneDir 	=  N * (2*(M+1)-1) - (N-1) ;
        nmbPoints_oneDir_allSubdomain 	= length * nmbPoints_oneDir - (length-1) ;
        nmbPoints 			= (2*(M+1)-1)*(2*(M+1)-1) ;
    }
    else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Wrong FE-Type, either P1 or P2.");
    }
    
    this->FEType_ = FEType;
    
    GO nmbPGlob_oneDir = N * (M+1) - (N-1);
    
    this->numElementsGlob_ = 2*(nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1) * bfs_multiplier;
    
    if (rank>=size) {
        M = -1; // keine Schleife wird ausgefuehrt
        nmbElements = 0;
        nmbPoints = 0;
    }
    else{
        nmbElements = 2*(M)*(M);
    }
    
    // P0 Mesh
    if (FEType == "P0") {
        if (verbose)
            cout << "-- H:"<<H << " h:" <<h << " --" << endl;

        if (verbose)
            cout << "-- Building P0 Points Repeated ... " << endl;
        
        this->pointsRep_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(2,0.0)));
        this->bcFlagRep_.reset(new std::vector<int> (nmbPoints,0));
        elementsVec = Teuchos::rcp(new std::vector<std::vector<int> >(nmbElements, std::vector<int>(3,-1)));
        
        Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);
        
        int whichSquareSet = (int)rank / nmbSubdomainsSquares;
        
        int offset_Squares_x = (int) (whichSquareSet+1) / 2;
        int offset_Squares_y = ((whichSquareSet+1) % 2);
        int counter = 0;
        int offset_x = ((rank - nmbSubdomainsSquares*whichSquareSet) % N);
        int offset_y = 0;
        
        if (((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))>=N) {
            offset_y = (int) ((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))/(N);
        }
        
        for (int s=0; s < M+1; s++) {
            for (int r=0; r < M+1; r++) {
                (*this->pointsRep_)[counter][0] = coorRec[0] + r*h + offset_x * H + offset_Squares_x * H * nmbSubdomainsSquares_OneDir;
                (*this->pointsRep_)[counter][1] = coorRec[1] + s*h + offset_y * H + offset_Squares_y * H * nmbSubdomainsSquares_OneDir;
                pointsRepGlobMapping[counter] = r + s*(nmbPoints_oneDir_allSubdomain - (1-offset_Squares_y)*(nmbPoints_oneDir-1))
                + offset_x*(M)
                + offset_y*((nmbPoints_oneDir_allSubdomain) - (1-offset_Squares_y)*(nmbPoints_oneDir-1)) *M
                + offset_Squares_x * M * nmbSubdomainsSquares_OneDir
                + offset_Squares_y * (nmbPoints_oneDir_allSubdomain * M * nmbSubdomainsSquares_OneDir
                                      - (nmbPoints_oneDir-1)*(nmbPoints_oneDir-1));
                //- M * nmbSubdomainsSquares_OneDir);
                if (offset_Squares_x>0 && offset_Squares_y==0 ) {
                    pointsRepGlobMapping[counter] -= nmbPoints_oneDir-1;
                }
                if (offset_Squares_x>0 && offset_Squares_y==0 && offset_y+1==nmbSubdomainsSquares_OneDir && s==M) {
                    pointsRepGlobMapping[counter] += nmbPoints_oneDir-1;
                }
                if ((*this->pointsRep_)[counter][0] > (coorRec[0]+length-eps) 	|| (*this->pointsRep_)[counter][0] < (coorRec[0]+eps) ||
                    (*this->pointsRep_)[counter][1] > (coorRec[1]+height-eps) 	|| (*this->pointsRep_)[counter][1] < (coorRec[1]+eps)) {
                    (*this->bcFlagRep_)[counter] = 1;
                }
           
                counter++;
            }
        }
        
        if (verbose) {
            cout << " done! --" << endl;
        }
        
        if (verbose) {
            cout << "-- Building P0 Repeated and Unique Map ... " << flush;
        }
        
        

        this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );

        this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );
        if (verbose) {
            cout << " done! --" << endl;
        }
        
        if (verbose) {
            cout << "-- Building P0 Unique Points ... " << flush;
        }
        
        this->pointsUni_.reset(new std::vector<std::vector<double> >(this->mapUnique_->getNodeNumElements(), std::vector<double>(2,0.0)));
        this->bcFlagUni_.reset(new std::vector<int> (this->mapUnique_->getNodeNumElements(),0));
        LO index;
        for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
            
            index = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement(i) );
            
            (*this->pointsUni_)[i][0] = (*this->pointsRep_)[index][0];
            (*this->pointsUni_)[i][1] = (*this->pointsRep_)[index][1];
            (*this->bcFlagUni_)[i] = (*this->bcFlagRep_)[index];
        }
        
        if (verbose) {
            cout << " done! --" << endl;
        }
        
        
        if (verbose) {
            cout << "-- Building P0 Elements ... " << flush;
        }
        
        counter = 0;
        for (int s=0; s < M; s++) {
            for (int r=0; r < M; r++) {
                (*elementsVec)[counter][0] = r+1 + (M+1) * s;
                (*elementsVec)[counter][1] = r + (M+1) * s ;
                (*elementsVec)[counter][2] = r+1 + (M+1) * (s+1);
                counter++;
                (*elementsVec)[counter][0] = r + (M+1) * (s+1);
                (*elementsVec)[counter][1] = r + (M+1) * (s);
                (*elementsVec)[counter][2] = r+1 + (M+1) * (s+1);
                counter++;
            }
        }
        
        if (verbose) {
            cout << " done! --" << endl;
        }
    }
    
    // P1 Mesh
    else if (FEType == "P1") {
        if (verbose) {
            cout << "-- H:"<<H << " h:" <<h << " --" << endl;
        }
        if (verbose) {
            cout << "-- Building P1 Points Repeated ... " << endl;
        }
        
        this->pointsRep_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(2,0.0)));
        this->bcFlagRep_.reset(new std::vector<int> (nmbPoints,0));
        elementsVec = Teuchos::rcp(new std::vector<std::vector<int> >(nmbElements,std::vector<int>(3,-1)));

        Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);
        
        int whichSquareSet = (int)rank / nmbSubdomainsSquares;
        
        int offset_Squares_x = (int) (whichSquareSet+1) / 2;
        int offset_Squares_y = ((whichSquareSet+1) % 2);
        int counter = 0;
        int offset_x = ((rank - nmbSubdomainsSquares*whichSquareSet) % N);
        int offset_y = 0;
                
        if (((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))>=N) {
            offset_y = (int) ((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))/(N);
        }
        
        for (int s=0; s < M+1; s++) {
            for (int r=0; r < M+1; r++) {
                (*this->pointsRep_)[counter][0] = coorRec[0] + r*h + offset_x * H + offset_Squares_x * H * nmbSubdomainsSquares_OneDir;
                (*this->pointsRep_)[counter][1] = coorRec[1] + s*h + offset_y * H + offset_Squares_y * H * nmbSubdomainsSquares_OneDir;
                pointsRepGlobMapping[counter] = r + s*(nmbPoints_oneDir_allSubdomain - (1-offset_Squares_y)*(nmbPoints_oneDir-1))
                + offset_x*(M)
                + offset_y*((nmbPoints_oneDir_allSubdomain) - (1-offset_Squares_y)*(nmbPoints_oneDir-1)) *M
                + offset_Squares_x * M * nmbSubdomainsSquares_OneDir
                + offset_Squares_y * (nmbPoints_oneDir_allSubdomain * M * nmbSubdomainsSquares_OneDir
                                      - (nmbPoints_oneDir-1)*(nmbPoints_oneDir-1));
                //- M * nmbSubdomainsSquares_OneDir);
                if (offset_Squares_x>0 && offset_Squares_y==0 ) {
                    pointsRepGlobMapping[counter] -= nmbPoints_oneDir-1;
                }
                if (offset_Squares_x>0 && offset_Squares_y==0 && offset_y+1==nmbSubdomainsSquares_OneDir && s==M) {
                    pointsRepGlobMapping[counter] += nmbPoints_oneDir-1;
                }
                if ((*this->pointsRep_)[counter][0] > (coorRec[0]+length-eps) 	|| (*this->pointsRep_)[counter][0] < (coorRec[0]+eps) ||
                    (*this->pointsRep_)[counter][1] > (coorRec[1]+height-eps) 	|| (*this->pointsRep_)[counter][1] < (coorRec[1]+eps)) {
                    (*this->bcFlagRep_)[counter] = 1;
                }
                
                counter++;
            }
        }
        
        if (verbose) {
            cout << " done! --" << endl;
        }
        
        if (verbose) {
            cout << "-- Building P1 Repeated and Unique Map ... " << flush;
        }
        
        

    
        this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );

        this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );
        
        if (verbose) {
            cout << " done! --" << endl;
        }
        
        if (verbose) {
            cout << "-- Building P1 Unique Points ... " << flush;
        }
        
        this->pointsUni_.reset(new std::vector<std::vector<double> >(this->mapUnique_->getNodeNumElements(), std::vector<double>(2,0.0)));
        this->bcFlagUni_.reset(new std::vector<int> (this->mapUnique_->getNodeNumElements(),0));
        LO index;
        for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
            
            index = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement(i) );
            
            (*this->pointsUni_)[i][0] = (*this->pointsRep_)[index][0];
            (*this->pointsUni_)[i][1] = (*this->pointsRep_)[index][1];
            (*this->bcFlagUni_)[i] = (*this->bcFlagRep_)[index];
        }
        
        if (verbose) {
            cout << " done! --" << endl;
        }
        
        
        if (verbose) {
            cout << "-- Building P1 Elements ... " << flush;
        }
        
        counter = 0;
        for (int s=0; s < M; s++) {
            for (int r=0; r < M; r++) {
                (*elementsVec)[counter][0] = r+1 + (M+1) * s;
                (*elementsVec)[counter][1] = r + (M+1) * s ;
                (*elementsVec)[counter][2] = r+1 + (M+1) * (s+1);
                counter++;
                (*elementsVec)[counter][0] = r + (M+1) * (s+1);
                (*elementsVec)[counter][1] = r + (M+1) * (s);
                (*elementsVec)[counter][2] = r+1 + (M+1) * (s+1);
                counter++;
            }
        }
        
        if (verbose) {
            cout << " done! --" << endl;
        }
    }
    
    // P2 Mesh
    else if(FEType == "P2"){
        if (verbose) {
            cout << "-- H:"<<H << " h:" <<h << " --" << endl;
        }
        if (verbose) {
            cout << "-- Building P2 Points Repeated ... " << flush;
        }
        
        this->pointsRep_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(2,0.0)));
        this->bcFlagRep_.reset(new std::vector<int> (nmbPoints,0));
        elementsVec = Teuchos::rcp(new std::vector<std::vector<int> >(nmbElements,std::vector<int>(6,-1)));
        Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);
        
        int whichSquareSet = (int)rank / nmbSubdomainsSquares;
        
        int offset_Squares_x = (int) (whichSquareSet+1) / 2;
        int offset_Squares_y = ((whichSquareSet+1) % 2);
        
        int counter = 0;
        int offset_x = ((rank - nmbSubdomainsSquares*whichSquareSet) % N);
        int offset_y = 0;
        
        
        if (((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))>=N) {
            offset_y = (int) ((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))/(N);
        }
        
        
        bool p1point;
        int p1_s = 0;
        int p1_r = 0;
        
        for (int s=0; s < 2*(M+1)-1; s++) {
            for (int r=0; r < 2*(M+1)-1; r++) {
                p1point = false;
                if (s%2==0 && r%2==0) {
                    p1point = true;
                    p1_s = s/2;
                    p1_r = r/2;
                }
                (*this->pointsRep_)[counter][0] = coorRec[0] + r*h/2. + offset_x * H + offset_Squares_x * H * nmbSubdomainsSquares_OneDir;
                (*this->pointsRep_)[counter][1] = coorRec[1] + s*h/2. + offset_y * H + offset_Squares_y * H * nmbSubdomainsSquares_OneDir;
                pointsRepGlobMapping[counter] = r + s*(nmbPoints_oneDir_allSubdomain - (1-offset_Squares_y)*(nmbPoints_oneDir-1))
                + offset_x*(2*(M+1)-2)
                + offset_y*((nmbPoints_oneDir_allSubdomain) - (1-offset_Squares_y)*(nmbPoints_oneDir-1)) * (2*(M+1)-2)
                + offset_Squares_x * (2*(M+1)-2) * nmbSubdomainsSquares_OneDir
                + offset_Squares_y * (nmbPoints_oneDir_allSubdomain * 2*M * nmbSubdomainsSquares_OneDir
                                      - (nmbPoints_oneDir-1)*(nmbPoints_oneDir-1));
                //- M * nmbSubdomainsSquares_OneDir);
                if (offset_Squares_x>0 && offset_Squares_y==0 ) {
                    pointsRepGlobMapping[counter] -= nmbPoints_oneDir-1;
                }
                if (offset_Squares_x>0 && offset_Squares_y==0 && offset_y+1==nmbSubdomainsSquares_OneDir && s==2*M) {
                    pointsRepGlobMapping[counter] += nmbPoints_oneDir-1;
                }
                if ((*this->pointsRep_)[counter][0] > (coorRec[0]+length-eps) 	|| (*this->pointsRep_)[counter][0] < (coorRec[0]+eps) ||
                    (*this->pointsRep_)[counter][1] > (coorRec[1]+height-eps) 	|| (*this->pointsRep_)[counter][1] < (coorRec[1]+eps)) {
                    (*this->bcFlagRep_)[counter] = 1;
                }
                
                counter++;
            }
        }
        
        if (verbose) {
            cout << " done! --" << endl;
        }
        //        int* globIndex = new int[MappingPointsRepGlob->size()];
        //        globIndex = &(MappingPointsRepGlob->at(0));
        
        if (verbose) {
            cout << "-- Building P2 Repeated and Unique Map ... " << flush;
        }
        

        this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );

        this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );

        
        if (verbose) {
            cout << " done! --" << endl;
        }
        
        if (verbose) {
            cout << "-- Building P2 Unique Points ... " << flush;
        }
        
        this->pointsUni_.reset(new std::vector<std::vector<double> >(this->mapUnique_->getNodeNumElements(),std::vector<double>(2,0.0)));
        this->bcFlagUni_.reset(new std::vector<int> (this->mapUnique_->getNodeNumElements(),0));
        
        LO index;
        for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
            
            index = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement(i) );
            
            (*this->pointsUni_)[i][0] = (*this->pointsRep_)[index][0];
            (*this->pointsUni_)[i][1] = (*this->pointsRep_)[index][1];
            (*this->bcFlagUni_)[i] = (*this->bcFlagRep_)[index];
        }
        
        if (verbose) {
            cout << " done! --" << endl;
        }
        
        //                Triangle numbering
        //                    2
        //                  * *
        //                *   *
        //              4	  5
        //            *       *
        //          *         *
        //        1 * * 3 * * 0
        
        
        if (verbose) {
            cout << "-- Building P2 Elements ... " << flush;
        }
        
        int    P2M = 2*(M+1)-1;
        
        counter = 0;
        
        for (int s=0; s < M; s++) {
            for (int r=0; r < M; r++) {
                
                (*elementsVec)[counter][0] = 2*(r+1)	+ 2*P2M * (s) ;
                (*elementsVec)[counter][1] = 2*(r) 	+ 2*P2M * (s) ;
                (*elementsVec)[counter][2] = 2*(r+1)	+ 2*P2M * (s+1) ;
                
                (*elementsVec)[counter][3] = 2*(r) +1	+ 2*P2M * (s) ;
                (*elementsVec)[counter][4] = 2*(r) +1	+ 2*P2M * (s) +P2M ;
                (*elementsVec)[counter][5] = 2*(r+1)		+ 2*P2M * (s) +P2M ;
                
                counter++;
                
                (*elementsVec)[counter][0] = 2*(r) 	+ 2*P2M * (s+1) ;
                (*elementsVec)[counter][1] = 2*(r) 	+ 2*P2M * (s) ;
                (*elementsVec)[counter][2] = 2*(r+1)	+ 2*P2M * (s+1) ;
                
                (*elementsVec)[counter][3] = 2*(r)		+ 2*P2M * (s) +P2M ;
                (*elementsVec)[counter][4] = 2*(r) +1 	+ 2*P2M * (s) +P2M ;
                (*elementsVec)[counter][5] = 2*(r) +1	+ 2*P2M * (s+1) ;
                
                counter++;
                
            }
        }
        
        if (verbose) {
            cout << " done! --" << endl;
        }
    }
    buildElementsClass(elementsVec);
}

template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::buildMesh3DBFS(std::string FEType,
                                                    int N,
                                                    int M,
                                                    int numProcsCoarseSolve,
                                                    std::string underlyingLib){
    
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ScalarTraits;
    
    TEUCHOS_TEST_FOR_EXCEPTION(!(M>=1),std::logic_error,"H/h is to small.");
    TEUCHOS_TEST_FOR_EXCEPTION(this->comm_.is_null(),std::runtime_error,"comm_ is null.");
    
    typedef ScalarTraits<SC> ST;
    SC eps = ST::eps();
    
    bool verbose (this->comm_->getRank() == 0);
    
    setRankRange( numProcsCoarseSolve );
    
    int         rank = this->comm_->getRank();
    int         size = this->comm_->getSize() - numProcsCoarseSolve;

    int         bfs_multiplier = (int) 2*(length)-1;

    int         nmbSubdomainsSquares = size / bfs_multiplier;
    int         nmbSubdomainsSquares_OneDir = (std::pow(nmbSubdomainsSquares,1./3.) + 100*eps); // same as N
    
    SC      h = ST::one()/(M*N);
    SC      H = ST::one()/N;
    
    LO nmbElements;
    LO nmbPoints;

    GO   nmbPoints_oneDir;
    GO   nmbPoints_oneDir_allSubdomain;

    if (FEType == "P0") {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"implement P0.");
    }
    else if (FEType == "P1") {
        nmbPoints_oneDir 	= N * (M+1) - (N-1);
        nmbPoints_oneDir_allSubdomain = length * nmbPoints_oneDir - (length-1);
        nmbPoints			= (M+1)*(M+1)*(M+1);
    }
    else if(FEType == "P2"){
        nmbPoints_oneDir 	=  N * (2*(M+1)-1) - (N-1);
        nmbPoints_oneDir_allSubdomain = length * nmbPoints_oneDir - (length-1) ;
        nmbPoints 			= (2*(M+1)-1)*(2*(M+1)-1)*(2*(M+1)-1);
    }
    else if(FEType == "P2-CR"){
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"P2-CR might not work properly.");        
        nmbPoints_oneDir 	=  N * (2*(M+1)-1) - (N-1);
        nmbPoints_oneDir_allSubdomain = length * nmbPoints_oneDir - (length-1) ;
        nmbPoints 			= (2*(M+1)-1)*(2*(M+1)-1)*(2*(M+1)-1);
    }
    else if(FEType == "P1-disc" || FEType == "P1-disc-global"){
        nmbPoints_oneDir 	= N * (M+1) - (N-1);
        nmbPoints_oneDir_allSubdomain = length * nmbPoints_oneDir - (length-1);
        nmbPoints			= (M+1)*(M+1)*(M+1);
    }
    else if(FEType == "Q2"){
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Wrong FE-Type, only P1,P1-disc, P1-disc-global, P2, or P2-CR.");
    
    this->FEType_ = FEType;
    
    GO nmbPGlob_oneDir = N * (M+1) - (N-1);
    this->numElementsGlob_ = 6*(nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1)*(nmbPGlob_oneDir-1) * bfs_multiplier;
    int MM=M;
    if (rank>=size) {
        M = -1; // keine Schleife wird ausgefuehrt
        nmbElements = 0;
        nmbPoints = 0;
    }
    else{
        nmbElements = 6*M*M*M;
    }
    
    // P1 Mesh
    if (FEType == "P1") {
        
        this->pointsRep_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(3,0.0)));
        this->bcFlagRep_.reset(new std::vector<int> (nmbPoints,0));
        vec2D_int_ptr_Type elementsVec = Teuchos::rcp(new std::vector<std::vector<int> >(nmbElements,std::vector<int>(4,-1)));

        Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);
        
        int whichSquareSet = (int)rank / nmbSubdomainsSquares;
        
        int offset_Squares_x = (int) (whichSquareSet+1) / 2;
        int offset_Squares_y = 0;
        int offset_Squares_z = ((whichSquareSet+1) % 2);
        
        int counter = 0;
        int offset_x = ((rank - nmbSubdomainsSquares*whichSquareSet) % N);
        int offset_y = 0;
        int offset_z = 0;
        if (((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))>=N)
            offset_y = (int) ((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))/(N);
        if (((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N*N))>=N*N )
            offset_z = (int) ((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N*N))/(N*(N));
        
        for (int t=0; t < M+1; t++) {
            for (int s=0; s < M+1; s++) {
                for (int r=0; r < M+1; r++) {
                    (*this->pointsRep_)[counter][0] = coorRec[0] + r*h + offset_x * H + offset_Squares_x * H * nmbSubdomainsSquares_OneDir;
                    
                    (*this->pointsRep_)[counter][1] = coorRec[1] + s*h + offset_y * H + offset_Squares_y * H * nmbSubdomainsSquares_OneDir;
                    
                    (*this->pointsRep_)[counter][2] = coorRec[2] + t*h + offset_z * H + offset_Squares_z * H * nmbSubdomainsSquares_OneDir;
                    
                    pointsRepGlobMapping[counter] = r
                    + s*(nmbPoints_oneDir_allSubdomain);
                    if (offset_Squares_x > 0 && offset_Squares_z == 0 && offset_z+1!=nmbSubdomainsSquares_OneDir) {
                        pointsRepGlobMapping[counter] -= s*(nmbPoints_oneDir-1);
                    }
                    else if(offset_Squares_x > 0 && offset_Squares_z == 0 && offset_z+1==nmbSubdomainsSquares_OneDir && t!=M){
                        pointsRepGlobMapping[counter] -= s*(nmbPoints_oneDir-1);
                    }
                    
                    pointsRepGlobMapping[counter] += t*nmbPoints_oneDir_allSubdomain*nmbPoints_oneDir;
                    if (offset_Squares_x > 0 && offset_Squares_z == 0 ) {
                        pointsRepGlobMapping[counter] -= t*(nmbPoints_oneDir-1)*(nmbPoints_oneDir);
                    }
                    
                    pointsRepGlobMapping[counter] += offset_x*(M)
                    + offset_y*( nmbPoints_oneDir_allSubdomain * M );
                    if (offset_Squares_x > 0 && offset_Squares_z == 0 && offset_z+1!=nmbSubdomainsSquares_OneDir) {
                        pointsRepGlobMapping[counter] -= offset_y*M*(nmbPoints_oneDir-1);
                    }
                    else if(offset_Squares_x > 0 && offset_Squares_z == 0 && offset_z+1==nmbSubdomainsSquares_OneDir && t!=M){
                        pointsRepGlobMapping[counter] -= offset_y*M*(nmbPoints_oneDir-1);
                    }
                    
                    pointsRepGlobMapping[counter] += offset_z * M * nmbPoints_oneDir_allSubdomain * nmbPoints_oneDir;
                    if (offset_Squares_x > 0 && offset_Squares_z == 0 ) {
                        pointsRepGlobMapping[counter] -= offset_z*M*(nmbPoints_oneDir-1)*(nmbPoints_oneDir);
                    }
                    
                    pointsRepGlobMapping[counter] += offset_Squares_x * M * nmbSubdomainsSquares_OneDir;
                    if (offset_Squares_z == 0 && offset_z+1!=nmbSubdomainsSquares_OneDir) {
                        pointsRepGlobMapping[counter] -= M * nmbSubdomainsSquares_OneDir;
                    }
                    else if(offset_Squares_z == 0 && offset_z+1==nmbSubdomainsSquares_OneDir && t!=M){
                        pointsRepGlobMapping[counter] -= M * nmbSubdomainsSquares_OneDir;
                    }
                    
                    pointsRepGlobMapping[counter] += offset_Squares_z * nmbPoints_oneDir_allSubdomain * ((M) * nmbSubdomainsSquares_OneDir+1) * M * nmbSubdomainsSquares_OneDir;
                    if (offset_Squares_z > 0 ) {
                        pointsRepGlobMapping[counter] -= (nmbPoints_oneDir-1)*(nmbPoints_oneDir-1)*(nmbPoints_oneDir);
                    }
                    counter++;
                }
            }
        }
        

        this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );
        
        this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );
        
        if (verbose) {
            cout << "-- Building P1 Unique Points ... " << flush;
        }
        
        this->pointsUni_.reset(new std::vector<std::vector<double> >(this->mapUnique_->getNodeNumElements(), std::vector<double>(3,0.0)));
        this->bcFlagUni_.reset(new std::vector<int> (this->mapUnique_->getNodeNumElements(),0));
        LO index;

        for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
            index = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement(i) );
            (*this->pointsUni_)[i][0] = (*this->pointsRep_)[index][0];
            (*this->pointsUni_)[i][1] = (*this->pointsRep_)[index][1];
            (*this->pointsUni_)[i][2] = (*this->pointsRep_)[index][2];
            (*this->bcFlagUni_)[i] = (*this->bcFlagRep_)[index];
        }

        if (verbose) {
            cout << " done! --" << endl;
        }
        
        
        if (verbose) {
            cout << "-- Building P1 Elements ... " << flush;
        }
        counter = 0;
        for (int t=0; t < M; t++) {
            for (int s=0; s < M; s++) {
                for (int r=0; r < M; r++) {
                    (*elementsVec)[counter][0] = r+1 + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][1] = r + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][2] = r+1 + (M+1) * s + (M+1)*(M+1) * (t+1) ;
                    (*elementsVec)[counter][3] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    counter++;
                    (*elementsVec)[counter][0] = r + (M+1) * s + (M+1)*(M+1) * (t+1) ;
                    (*elementsVec)[counter][1] = r + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][2] = r+1 + (M+1) * s + (M+1)*(M+1) * (t+1) ;
                    (*elementsVec)[counter][3] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    counter++;
                    (*elementsVec)[counter][0] = r+1 + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][1] = r + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][2] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][3] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    counter++;
                    (*elementsVec)[counter][0] = r + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][1] = r + (M+1) * (s+1) + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][2] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][3] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    counter++;
                    (*elementsVec)[counter][0] = r + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][1] = r + (M+1) * (s+1) + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][2] = r + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    (*elementsVec)[counter][3] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    counter++;
                    (*elementsVec)[counter][0] = r + (M+1) * s + (M+1)*(M+1) * t ;
                    (*elementsVec)[counter][1] = r + (M+1) * s + (M+1)*(M+1) * (t+1) ;
                    (*elementsVec)[counter][2] = r + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    (*elementsVec)[counter][3] = r+1 + (M+1) * (s+1) + (M+1)*(M+1) * (t+1) ;
                    counter++;
                }
            }
        }
        if (verbose) {
            cout << " done! --" << endl;
        }
        buildElementsClass(elementsVec);
    }
    else if(FEType == "P1-disc" || FEType == "P1-disc-global")
        buildP1_Disc_Q2_3DBFS( N, MM, numProcsCoarseSolve, underlyingLib );
    else if(FEType == "P2"){
        
        this->pointsRep_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(3,0.0)));
        this->bcFlagRep_.reset(new std::vector<int> (nmbPoints,0));
        vec2D_int_ptr_Type elementsVec = Teuchos::rcp(new std::vector<std::vector<int> >(nmbElements,std::vector<int>(10,-1)));
        Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);
        
        int whichSquareSet = (int)rank / nmbSubdomainsSquares;
        
        int offset_Squares_x = (int) (whichSquareSet+1) / 2;
        int offset_Squares_y = 0;
        int offset_Squares_z = ((whichSquareSet+1) % 2);
        
        int counter = 0;
        int offset_x = ((rank - nmbSubdomainsSquares*whichSquareSet) % N);
        int offset_y = 0;
        int offset_z = 0;
        if (((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))>=N) {
            offset_y = (int) ((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))/(N);
        }
        
        if (((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N*N))>=N*N ) {
            offset_z = (int) ((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N*N))/(N*(N));
        }
        
        bool p1point;
        int p1_s = 0;
        int p1_r = 0;
        int p1_t = 0;
        for (int t=0; t < 2*(M+1)-1; t++) {
            for (int s=0; s < 2*(M+1)-1; s++) {
                for (int r=0; r < 2*(M+1)-1; r++) {
                    p1point = false;
                    if (s%2==0 && r%2==0 && t%2==0) {
                        p1point = true;
                        p1_s = s/2;
                        p1_r = r/2;
                        p1_t = t/2;
                    }
                    (*this->pointsRep_)[counter][0] = coorRec[0] + r*h/2. + offset_x * H + offset_Squares_x * H * nmbSubdomainsSquares_OneDir;
                    
                    (*this->pointsRep_)[counter][1] = coorRec[1] + s*h/2. + offset_y * H + offset_Squares_y * H * nmbSubdomainsSquares_OneDir;
                    
                    (*this->pointsRep_)[counter][2] = coorRec[2] + t*h/2. + offset_z * H + offset_Squares_z * H * nmbSubdomainsSquares_OneDir;
                    
                    pointsRepGlobMapping[counter] = r
                    + s*(nmbPoints_oneDir_allSubdomain);
                    if (offset_Squares_x > 0 && offset_Squares_z == 0 && offset_z+1!=nmbSubdomainsSquares_OneDir) {
                        pointsRepGlobMapping[counter] -= s*(nmbPoints_oneDir-1);
                    }
                    else if(offset_Squares_x > 0 && offset_Squares_z == 0 && offset_z+1==nmbSubdomainsSquares_OneDir && t!=2*M){
                        pointsRepGlobMapping[counter] -= s*(nmbPoints_oneDir-1);
                    }
                    
                    pointsRepGlobMapping[counter] += t*nmbPoints_oneDir_allSubdomain*nmbPoints_oneDir;
                    if (offset_Squares_x > 0 && offset_Squares_z == 0 ) {
                        pointsRepGlobMapping[counter] -= t*(nmbPoints_oneDir-1)*(nmbPoints_oneDir);
                    }
                    
                    pointsRepGlobMapping[counter] += offset_x*(2*M)
                    + offset_y*( nmbPoints_oneDir_allSubdomain * 2*M );
                    if (offset_Squares_x > 0 && offset_Squares_z == 0 && offset_z+1!=nmbSubdomainsSquares_OneDir) {
                        pointsRepGlobMapping[counter] -= offset_y*2*M*(nmbPoints_oneDir-1);
                    }
                    else if(offset_Squares_x > 0 && offset_Squares_z == 0 && offset_z+1==nmbSubdomainsSquares_OneDir && t!=2*M){
                        pointsRepGlobMapping[counter] -= offset_y*2*M*(nmbPoints_oneDir-1);
                    }
                    
                    pointsRepGlobMapping[counter] += offset_z * 2*M * nmbPoints_oneDir_allSubdomain * nmbPoints_oneDir;
                    if (offset_Squares_x > 0 && offset_Squares_z == 0 ) {
                        pointsRepGlobMapping[counter] -= offset_z*2*M*(nmbPoints_oneDir-1)*(nmbPoints_oneDir);
                    }
                    
                    pointsRepGlobMapping[counter] += offset_Squares_x * 2*M * nmbSubdomainsSquares_OneDir;
                    if (offset_Squares_z == 0 && offset_z+1!=nmbSubdomainsSquares_OneDir) {
                        pointsRepGlobMapping[counter] -= 2*M * nmbSubdomainsSquares_OneDir;
                    }
                    else if(offset_Squares_z == 0 && offset_z+1==nmbSubdomainsSquares_OneDir && t!=2*M){
                        pointsRepGlobMapping[counter] -= 2*M * nmbSubdomainsSquares_OneDir;
                    }
                    
                    pointsRepGlobMapping[counter] += offset_Squares_z * nmbPoints_oneDir_allSubdomain * ((2*M) * nmbSubdomainsSquares_OneDir+1) * 2*M * nmbSubdomainsSquares_OneDir;
                    if (offset_Squares_z > 0 ) {
                        pointsRepGlobMapping[counter] -= (nmbPoints_oneDir-1)*(nmbPoints_oneDir-1)*(nmbPoints_oneDir);
                    }
                    counter++;
                }
            }
        }
        
        this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );
        
        this->mapUnique_ = this->mapRepeated_->buildUniqueMap( numProcsCoarseSolve );
        
        if (verbose) {
            cout << "-- Building P2 Unique Points ... " << flush;
        }
        
        this->pointsUni_.reset(new std::vector<std::vector<double> >(this->mapUnique_->getNodeNumElements(),std::vector<double>(3,0.0)));
        this->bcFlagUni_.reset(new std::vector<int> (this->mapUnique_->getNodeNumElements(),0));
        
        LO index;
        for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
            
            index = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement(i) );
            
            (*this->pointsUni_)[i][0] = (*this->pointsRep_)[index][0];
            (*this->pointsUni_)[i][1] = (*this->pointsRep_)[index][1];
            (*this->pointsUni_)[i][2] = (*this->pointsRep_)[index][2];
            (*this->bcFlagUni_)[i] = (*this->bcFlagRep_)[index];
        }

        if (verbose) {
            cout << " done! --" << endl;
        }
        
        //                Face 1          Face2               Face 3            Face 4
        //                    2      2 * * 9 * * 3        3 * * 9 * * 2          	3
        //                  * *      *          *          *          * 		  * *
        //                *   *      *        *             *        *          *   *
        //              5	  6      6      7                8      5         8	    7
        //            *       *      *    *                   *    *        *       *
        //          *         *      *  *                      *  *       *         *
        //        1 * * 4 * * 0       0                         1       1 * * 4 * * 0
        
        
        int    P2M = 2*(M+1)-1;
        
        counter = 0;
        for (int t=0; t < M; t++) {
            for (int s=0; s < M; s++) {
                for (int r=0; r < M; r++) {
                    
                    (*elementsVec)[counter][0] = 2*(r+1)	+ 2*P2M * (s) 	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][1] = 2*(r) 	+ 2*P2M * (s) 	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][2] = 2*(r+1)	+ 2*P2M * (s) 	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][3] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    
                    (*elementsVec)[counter][4] = 2*(r) +1	+ 2*P2M * (s)		+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][6] = 2*(r+1)		+ 2*P2M * (s)		+ 2*P2M*P2M * (t) +P2M*P2M;
                    (*elementsVec)[counter][7] = 2*(r+1)		+ 2*P2M * (s) +P2M 	+ 2*P2M*P2M * (t) +P2M*P2M;
                    (*elementsVec)[counter][5] = 2*(r) +1	+ 2*P2M * (s) 		+ 2*P2M*P2M * (t) +P2M*P2M;
                    (*elementsVec)[counter][8] = 2*(r) +1	+ 2*P2M * (s) +P2M 	+ 2*P2M*P2M * (t) +P2M*P2M;
                    (*elementsVec)[counter][9] = 2*(r+1)		+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t+1);
                    
                    counter++;
                    
                    (*elementsVec)[counter][0] = 2*(r) 	+ 2*P2M * (s)	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][1] = 2*(r) 	+ 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][2] = 2*(r+1)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][3] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    
                    (*elementsVec)[counter][4] = 2*(r) 		+ 2*P2M * (s) 		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][6] = 2*(r) +1	+ 2*P2M * (s) 		+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][7] = 2*(r) +1	+ 2*P2M * (s)+P2M	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][5] = 2*(r) +1	+ 2*P2M * (s) 		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][8] = 2*(r) +1	+ 2*P2M * (s)+P2M 	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][9] = 2*(r+1)		+ 2*P2M * (s)+P2M 	+ 2*P2M*P2M * (t+1);
                    
                    counter++;
                    
                    (*elementsVec)[counter][0] = 2*(r+1)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][1] = 2*(r)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][2] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][3] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    
                    (*elementsVec)[counter][4] = 2*(r) +1	+ 2*P2M * (s)		+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][6] = 2*(r+1)		+ 2*P2M * (s)+P2M	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][7] = 2*(r+1)		+ 2*P2M * (s)+P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][5] = 2*(r) +1	+ 2*P2M * (s)+P2M	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][8] = 2*(r) +1	+ 2*P2M * (s)+P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][9] = 2*(r+1)		+ 2*P2M * (s+1)		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    
                    counter++;
                    
                    (*elementsVec)[counter][0] = 2*(r)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][1] = 2*(r)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][2] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][3] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    
                    (*elementsVec)[counter][4] = 2*(r) 		+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][6] = 2*(r) +1	+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][7] = 2*(r) +1	+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][5] = 2*(r) +1	+ 2*P2M * (s+1)		+ 2*P2M*P2M * (t);
                    (*elementsVec)[counter][8] = 2*(r) +1	+ 2*P2M * (s+1)		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][9] = 2*(r+1)		+ 2*P2M * (s+1)		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    
                    counter++;
                    
                    (*elementsVec)[counter][0] = 2*(r)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][1] = 2*(r)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][2] = 2*(r)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][3] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    
                    (*elementsVec)[counter][4] = 2*(r)		+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][6] = 2*(r)		+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][7] = 2*(r) +1	+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][5] = 2*(r)		+ 2*P2M*(s+1) 		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][8] = 2*(r) +1	+ 2*P2M*(s+1) 		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][9] = 2*(r) +1 	+ 2*P2M*(s+1) 		+ 2*P2M*P2M * (t+1) ;
                    
                    counter++;
                    
                    (*elementsVec)[counter][0] = 2*(r)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t) ;
                    (*elementsVec)[counter][1] = 2*(r)	+ 2*P2M * (s)	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][2] = 2*(r)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][3] = 2*(r+1)	+ 2*P2M * (s+1)	+ 2*P2M*P2M * (t+1) ;
                    
                    (*elementsVec)[counter][4] = 2*(r)		+ 2*P2M * (s)		+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][6] = 2*(r)		+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][7] = 2*(r) +1	+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t) +P2M*P2M ;
                    (*elementsVec)[counter][5] = 2*(r) 		+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][8] = 2*(r) +1 	+ 2*P2M * (s) +P2M	+ 2*P2M*P2M * (t+1) ;
                    (*elementsVec)[counter][9] = 2*(r) +1	+ 2*P2M*(s+1)		+ 2*P2M*P2M * (t+1) ;
                    
                    counter++;
                    
                }
            }
        }
        buildElementsClass(elementsVec);
    }
    else if(FEType == "Q2")
        build3DQ2BFS( N, MM, numProcsCoarseSolve, underlyingLib );
    
};
    
template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::buildP1_Disc_Q2_3DBFS(int N,
                                                     int M,
                                                     int numProcsCoarseSolve,
                                                     std::string underlyingLib){

    
    
    
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ScalarTraits;
    typedef ScalarTraits<SC> ST;
    
    bool verbose (this->comm_->getRank() == 0);
    
    setRankRange( numProcsCoarseSolve );
    
    if (verbose)
        cout << endl;
    
    SC eps = ScalarTraits<SC>::eps();
    
    int         rank = this->comm_->getRank();
    int         size = this->comm_->getSize() - numProcsCoarseSolve;
    
    
    int         bfs_multiplier = (int) 2*(length)-1;
    
    int         nmbSubdomainsSquares = size / bfs_multiplier;
    int         nmbSubdomainsSquares_OneDir = (std::pow(nmbSubdomainsSquares,1./3.) + 100*eps); // same as N
    
    SC      h = ST::one()/(M*N);
    SC      H = ST::one()/N;
    
    LO nmbElements = M*M*M;
    LO nmbPoints = 4*M*M*M; // 4 points for each element
    

    if (rank>=size) {
        M = -1; // keine Schleife wird ausgefuehrt
        nmbElements = 0;
        nmbPoints = 0;
    }
    
    this->numElementsGlob_ = nmbElements * size;
    
    int whichSquareSet = (int)rank / nmbSubdomainsSquares;
    
    int offset_Squares_x = (int) (whichSquareSet+1) / 2;
    int offset_Squares_y = 0;
    int offset_Squares_z = ((whichSquareSet+1) % 2);
    
    int counter = 0;
    int offset_x = ((rank - nmbSubdomainsSquares*whichSquareSet) % N);
    int offset_y = 0;
    int offset_z = 0;
    if (((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))>=N)
        offset_y = (int) ((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N))/(N);
    
    if (((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N*N))>=N*N )
        offset_z = (int) ((rank - nmbSubdomainsSquares*whichSquareSet) % (N*N*N))/(N*(N));
    
    
    this->pointsRep_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(3,0.0)));
    this->pointsUni_.reset(new std::vector<std::vector<double> >(nmbPoints,std::vector<double>(3,0.0)));
    this->bcFlagRep_.reset(new std::vector<int> (nmbPoints,0));
    this->bcFlagUni_.reset(new std::vector<int> (nmbPoints,0));
    Teuchos::Array<GO> pointsRepGlobMapping(nmbPoints);           
    
    if (verbose)
        cout << "-- Building P1-disc Points and Elements according to Q2 ... " << flush;
    
    vec2D_int_ptr_Type elementsVec = Teuchos::rcp(new std::vector<std::vector<int> >(nmbElements,std::vector<int>(4,-1)));
    
    counter = 0;
    LO pCounter = 0;
    GO globalCounterPoints = rank * nmbPoints;
    for (int t=0; t < M; t++) {
        for (int s=0; s < M; s++) {
            for (int r=0; r < M; r++) {
                
                //point 1
                (*this->pointsRep_)[pCounter][0] = coorRec[0] + r*h + offset_x * H + offset_Squares_x * H * nmbSubdomainsSquares_OneDir;
                (*this->pointsRep_)[pCounter][1] = coorRec[1] + s*h + offset_y * H + offset_Squares_y * H * nmbSubdomainsSquares_OneDir;
                (*this->pointsRep_)[pCounter][2] = coorRec[2] + t*h + offset_z * H + offset_Squares_z * H * nmbSubdomainsSquares_OneDir;
                
                if ((*this->pointsRep_)[pCounter][0] > (coorRec[0]+length-eps) 	|| (*this->pointsRep_)[pCounter][0] < (coorRec[0]+eps) ||
                    (*this->pointsRep_)[pCounter][1] > (coorRec[1]+width-eps) 	|| (*this->pointsRep_)[pCounter][1] < (coorRec[1]+eps) ||
                    (*this->pointsRep_)[pCounter][2] > (coorRec[2]+height-eps) 	|| (*this->pointsRep_)[pCounter][2] < (coorRec[2]+eps) ) {
                    
                    (*this->bcFlagRep_)[counter] = 1;
                }
                
                (*this->pointsUni_)[pCounter][0] = (*this->pointsRep_)[pCounter][0];
                (*this->pointsUni_)[pCounter][1] = (*this->pointsRep_)[pCounter][1];
                (*this->pointsUni_)[pCounter][2] = (*this->pointsRep_)[pCounter][2];
                
                (*this->bcFlagUni_)[pCounter] = (*this->bcFlagRep_)[pCounter];
                
                (*elementsVec)[counter][0] = pCounter;
                pointsRepGlobMapping[pCounter] = globalCounterPoints;
                globalCounterPoints++;
                pCounter++;
                
                //point 2
                (*this->pointsRep_)[pCounter][0] = coorRec[0] + (r+1)*h + offset_x * H + offset_Squares_x * H * nmbSubdomainsSquares_OneDir;
                (*this->pointsRep_)[pCounter][1] = coorRec[1] + s*h + offset_y * H + offset_Squares_y * H * nmbSubdomainsSquares_OneDir;
                (*this->pointsRep_)[pCounter][2] = coorRec[2] + t*h + offset_z * H + offset_Squares_z * H * nmbSubdomainsSquares_OneDir;
                if ((*this->pointsRep_)[pCounter][0] > (coorRec[0]+length-eps) 	|| (*this->pointsRep_)[pCounter][0] < (coorRec[0]+eps) ||
                    (*this->pointsRep_)[pCounter][1] > (coorRec[1]+width-eps) 	|| (*this->pointsRep_)[pCounter][1] < (coorRec[1]+eps) ||
                    (*this->pointsRep_)[pCounter][2] > (coorRec[2]+height-eps) 	|| (*this->pointsRep_)[pCounter][2] < (coorRec[2]+eps) ) {
                    
                    (*this->bcFlagRep_)[counter] = 1;
                }
                
                (*this->pointsUni_)[pCounter][0] = (*this->pointsRep_)[pCounter][0];
                (*this->pointsUni_)[pCounter][1] = (*this->pointsRep_)[pCounter][1];
                (*this->pointsUni_)[pCounter][2] = (*this->pointsRep_)[pCounter][2];
                
                (*this->bcFlagUni_)[pCounter] = (*this->bcFlagRep_)[pCounter];
                
                (*elementsVec)[counter][1] = pCounter;
                pointsRepGlobMapping[pCounter] = globalCounterPoints;
                globalCounterPoints++;
                pCounter++;
                //point 3
                (*this->pointsRep_)[pCounter][0] = coorRec[0] + r*h + offset_x * H + offset_Squares_x * H * nmbSubdomainsSquares_OneDir;
                (*this->pointsRep_)[pCounter][1] = coorRec[1] + (s+1)*h + offset_y * H + offset_Squares_y * H * nmbSubdomainsSquares_OneDir;
                (*this->pointsRep_)[pCounter][2] = coorRec[2] + t*h + offset_z * H + offset_Squares_z * H * nmbSubdomainsSquares_OneDir;
                if ((*this->pointsRep_)[pCounter][0] > (coorRec[0]+length-eps) 	|| (*this->pointsRep_)[pCounter][0] < (coorRec[0]+eps) ||
                    (*this->pointsRep_)[pCounter][1] > (coorRec[1]+width-eps) 	|| (*this->pointsRep_)[pCounter][1] < (coorRec[1]+eps) ||
                    (*this->pointsRep_)[pCounter][2] > (coorRec[2]+height-eps) 	|| (*this->pointsRep_)[pCounter][2] < (coorRec[2]+eps) ) {
                    
                    (*this->bcFlagRep_)[counter] = 1;
                }
                
                (*this->pointsUni_)[pCounter][0] = (*this->pointsRep_)[pCounter][0];
                (*this->pointsUni_)[pCounter][1] = (*this->pointsRep_)[pCounter][1];
                (*this->pointsUni_)[pCounter][2] = (*this->pointsRep_)[pCounter][2];
                
                (*this->bcFlagUni_)[pCounter] = (*this->bcFlagRep_)[pCounter];
                
                (*elementsVec)[counter][2] = pCounter;
                pointsRepGlobMapping[pCounter] = globalCounterPoints;
                globalCounterPoints++;
                pCounter++;
                
                //point 4
                (*this->pointsRep_)[pCounter][0] = coorRec[0] + r*h + offset_x * H + offset_Squares_x * H * nmbSubdomainsSquares_OneDir;
                (*this->pointsRep_)[pCounter][1] = coorRec[1] + s*h + offset_y * H + offset_Squares_y * H * nmbSubdomainsSquares_OneDir;
                (*this->pointsRep_)[pCounter][2] = coorRec[2] + (t+1)*h + offset_z * H + offset_Squares_z * H * nmbSubdomainsSquares_OneDir;
                if ((*this->pointsRep_)[pCounter][0] > (coorRec[0]+length-eps) 	|| (*this->pointsRep_)[pCounter][0] < (coorRec[0]+eps) ||
                    (*this->pointsRep_)[pCounter][1] > (coorRec[1]+width-eps) 	|| (*this->pointsRep_)[pCounter][1] < (coorRec[1]+eps) ||
                    (*this->pointsRep_)[pCounter][2] > (coorRec[2]+height-eps) 	|| (*this->pointsRep_)[pCounter][2] < (coorRec[2]+eps) ) {
                    
                    (*this->bcFlagRep_)[counter] = 1;
                }
                
                (*this->pointsUni_)[pCounter][0] = (*this->pointsRep_)[pCounter][0];
                (*this->pointsUni_)[pCounter][1] = (*this->pointsRep_)[pCounter][1];
                (*this->pointsUni_)[pCounter][2] = (*this->pointsRep_)[pCounter][2];
                
                (*this->bcFlagUni_)[pCounter] = (*this->bcFlagRep_)[pCounter];
                
                (*elementsVec)[counter][3] = pCounter;
                pointsRepGlobMapping[pCounter] = globalCounterPoints;
                globalCounterPoints++;
                pCounter++;
                
                counter++;
            }
        }
    }
    this->mapRepeated_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );
    
    this->mapUnique_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, pointsRepGlobMapping(), 0, this->comm_) );
    
    
    if (verbose)
        cout << "done!" << endl;
    
}
template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::setStructuredMeshFlags(int flagsOption,string FEType){
    
    double tol=1.e-12;
    
    switch (this->dim_) {
        case 2:
            switch (flagsOption) {
                case 0:
                    break;
                case 1: //Rectangle left inflow
                    for (int i=0; i<this->pointsUni_->size(); i++) {
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0] - tol) && this->pointsUni_->at(i).at(1) < (coorRec[1] + tol) ) {
                            this->bcFlagUni_->at(i) = 1;
                        }
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0] - tol) && this->pointsUni_->at(i).at(1) > (coorRec[1] + height - tol) ) {
                            this->bcFlagUni_->at(i) = 1;
                        }
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0]+length - tol) && this->pointsUni_->at(i).at(1) > (coorRec[1] + tol) && this->pointsUni_->at(i).at(1) < (coorRec[1] + height - tol)) {
                            this->bcFlagUni_->at(i) = 3; //outflow
                        }
                        if (this->pointsUni_->at(i).at(0) < (coorRec[0] +tol)) {
                            this->bcFlagUni_->at(i) = 2; //inflow
                        }
                    }
                    for (int i=0; i<this->pointsRep_->size(); i++) {
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0] - tol) && this->pointsRep_->at(i).at(1) < (coorRec[1] + tol) ) {
                            this->bcFlagRep_->at(i) = 1;
                        }
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0] - tol) && this->pointsRep_->at(i).at(1) > (coorRec[1] + height - tol) ) {
                            this->bcFlagRep_->at(i) = 1;
                        }
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0]+length - tol) && this->pointsRep_->at(i).at(1) > (coorRec[1] + tol) && this->pointsRep_->at(i).at(1) < (coorRec[1] + height - tol)) {
                            this->bcFlagRep_->at(i) = 3;
                        }
                        if (this->pointsRep_->at(i).at(0) < (coorRec[0] +tol)) {
                            this->bcFlagRep_->at(i) = 2;
                        }
                    }
                    break;
                case 2: //BFS
                    for (int i=0; i<this->pointsUni_->size(); i++) {
                        if (this->pointsUni_->at(i).at(0) < (coorRec[0] +tol) && this->pointsUni_->at(i).at(1) > (coorRec[1]+1. -tol) && this->pointsUni_->at(i).at(1) < (coorRec[1]+ height +tol)) {
                            this->bcFlagUni_->at(i) = 2;
                        }
                        if (this->pointsUni_->at(i).at(1) < (coorRec[1] + tol) || this->pointsUni_->at(i).at(1) > (coorRec[1]+height - tol) || this->pointsUni_->at(i).at(0) > (coorRec[0]+length - tol)) {
                            this->bcFlagUni_->at(i) = 1;
                        }
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0] - tol) && this->pointsUni_->at(i).at(0) < (coorRec[0]+1. + tol) && this->pointsUni_->at(i).at(1) > (coorRec[1] + 1. - tol) && this->pointsUni_->at(i).at(1) < (coorRec[1] + 1. + tol)) {
                            this->bcFlagUni_->at(i) = 1;
                        }
                        if (this->pointsUni_->at(i).at(1) > (coorRec[1] - tol) && this->pointsUni_->at(i).at(1) < (coorRec[1] + 1. + tol) && this->pointsUni_->at(i).at(0) > (coorRec[0] + 1. - tol) && this->pointsUni_->at(i).at(0) < (coorRec[0] + 1. + tol)) {
                            this->bcFlagUni_->at(i) = 1;
                        }
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0]+length - tol) && this->pointsUni_->at(i).at(1) > (coorRec[1] + tol) && this->pointsUni_->at(i).at(1) < (coorRec[1] + height - tol)) {
                            this->bcFlagUni_->at(i) = 3;
                        }
                    }
                    for (int i=0; i<this->pointsRep_->size(); i++) {
                        if (this->pointsRep_->at(i).at(0) < (coorRec[0] +tol) && this->pointsRep_->at(i).at(1) > (coorRec[1]+1. -tol) && this->pointsRep_->at(i).at(1) < (coorRec[1]+ height +tol)) {
                            this->bcFlagRep_->at(i) = 2;
                        }
                        if (this->pointsRep_->at(i).at(1) < (coorRec[1] + tol) || this->pointsRep_->at(i).at(1) > (coorRec[1]+height - tol) || this->pointsRep_->at(i).at(0) > (coorRec[0]+length - tol)) {
                            this->bcFlagRep_->at(i) = 1;
                        }
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0] - tol) && this->pointsRep_->at(i).at(0) < (coorRec[0]+1. + tol) && this->pointsRep_->at(i).at(1) > (coorRec[1] + 1. - tol) && this->pointsRep_->at(i).at(1) < (coorRec[1] + 1. + tol)) {
                            this->bcFlagRep_->at(i) = 1;
                        }
                        if (this->pointsRep_->at(i).at(1) > (coorRec[1] - tol) && this->pointsRep_->at(i).at(1) < (coorRec[1] + 1. + tol) && this->pointsRep_->at(i).at(0) > (coorRec[0] + 1. - tol) && this->pointsRep_->at(i).at(0) < (coorRec[0] + 1. + tol)) {
                            this->bcFlagRep_->at(i) = 1;
                        }
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0]+length - tol) && this->pointsRep_->at(i).at(1) > (coorRec[1] + tol) && this->pointsRep_->at(i).at(1) < (coorRec[1] + height - tol)) {
                            this->bcFlagRep_->at(i) = 3;
                        }
                    }
                    break;
                case 3: //tpm 
                    for (int i=0; i<this->pointsUni_->size(); i++) {
                        if (this->pointsUni_->at(i).at(0) < (coorRec[0] +tol)) {
                            this->bcFlagUni_->at(i) = 2; //left
                        }
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0] - tol) && this->pointsUni_->at(i).at(1) < (coorRec[1] + tol) ) {
                            this->bcFlagUni_->at(i) = 1; //bottom
                        }
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0] - tol) && this->pointsUni_->at(i).at(1) > (coorRec[1] + height - tol) ) {
                            this->bcFlagUni_->at(i) = 4; //top
                        }
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0]+length - tol) && this->pointsUni_->at(i).at(1) > (coorRec[1] + tol) && this->pointsUni_->at(i).at(1) < (coorRec[1] + height - tol)) {
                            this->bcFlagUni_->at(i) = 3; //right
                        }
                    }
                    for (int i=0; i<this->pointsRep_->size(); i++) {
                        if (this->pointsRep_->at(i).at(0) < (coorRec[0] +tol)) {
                            this->bcFlagRep_->at(i) = 2;
                        }
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0] - tol) && this->pointsRep_->at(i).at(1) < (coorRec[1] + tol) ) {
                            this->bcFlagRep_->at(i) = 1;
                        }
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0] - tol) && this->pointsRep_->at(i).at(1) > (coorRec[1] + height - tol) ) {
                            this->bcFlagRep_->at(i) = 4;
                        }
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0]+length - tol) && this->pointsRep_->at(i).at(1) > (coorRec[1] + tol) && this->pointsRep_->at(i).at(1) < (coorRec[1] + height - tol)) {
                            this->bcFlagRep_->at(i) = 3;
                        }
                    }
                    break;
                case 4: //mini tpm, we set all values manually
                    if (FEType == "P2") {
                        this->bcFlagUni_->at(0) = 1;
                        this->bcFlagUni_->at(1) = 2;
                        this->bcFlagUni_->at(2) = 2;
                        this->bcFlagUni_->at(3) = 2;
                        this->bcFlagUni_->at(4) = 4; // Dirichlet and lineload
                        this->bcFlagUni_->at(5) = 3;
                        this->bcFlagUni_->at(6) = 0;
                        this->bcFlagUni_->at(7) = 0;
                        this->bcFlagUni_->at(8) = 0;
                        this->bcFlagUni_->at(9) = 5; // lineload
                        this->bcFlagUni_->at(10) = 1;
                        this->bcFlagUni_->at(11) = 2;
                        this->bcFlagUni_->at(12) = 2;
                        this->bcFlagUni_->at(13) = 2;
                        this->bcFlagUni_->at(14) = 4;
                        
                        this->bcFlagRep_->at(0) = 1;
                        this->bcFlagRep_->at(1) = 2;
                        this->bcFlagRep_->at(2) = 2;
                        this->bcFlagRep_->at(3) = 2;
                        this->bcFlagRep_->at(4) = 4;
                        this->bcFlagRep_->at(5) = 3;
                        this->bcFlagRep_->at(6) = 0;
                        this->bcFlagRep_->at(7) = 0;
                        this->bcFlagRep_->at(8) = 0;
                        this->bcFlagRep_->at(9) = 5;
                        this->bcFlagRep_->at(10) = 1;
                        this->bcFlagRep_->at(11) = 2;
                        this->bcFlagRep_->at(12) = 2;
                        this->bcFlagRep_->at(13) = 2;
                        this->bcFlagRep_->at(14) = 4;
                    } else if(FEType=="P1") {
                        this->bcFlagUni_->at(0) = 1;
                        this->bcFlagUni_->at(1) = 1;
                        this->bcFlagUni_->at(2) = 1;
                        this->bcFlagUni_->at(3) = 2;
                        this->bcFlagUni_->at(4) = 2;
                        this->bcFlagUni_->at(5) = 1;
                        
                        this->bcFlagRep_->at(0) = 1;
                        this->bcFlagRep_->at(1) = 1;
                        this->bcFlagRep_->at(2) = 1;
                        this->bcFlagRep_->at(3) = 2;
                        this->bcFlagRep_->at(4) = 2;
                        this->bcFlagRep_->at(5) = 1;
                    }
                    break;
                default:
                    break;
            }
            break;
        case 3:
            switch (flagsOption) {
                case 0:
                    break;
                case 1: //Rectangle left inflow
                    for (int i=0; i<this->pointsUni_->size(); i++) {
                        if (this->pointsUni_->at(i).at(0) < (coorRec[0] + tol) ) {
                            this->bcFlagUni_->at(i) = 2;
                        }
                        //bottom
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0] + tol) &&
                            this->pointsUni_->at(i).at(2) < (coorRec[2] + tol) ) {
                            this->bcFlagUni_->at(i) = 1;
                        }
                        //top
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0] + tol) &&
                            this->pointsUni_->at(i).at(2) > (coorRec[2] + height - tol) ) {
                            this->bcFlagUni_->at(i) = 1;
                        }
                        //front
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0] + tol) &&
                            this->pointsUni_->at(i).at(1) < (coorRec[1] + tol) ) {
                            this->bcFlagUni_->at(i) = 1;
                        }
                        //back
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0] + tol) &&
                            this->pointsUni_->at(i).at(1) > (coorRec[1] + width - tol) ) {
                            this->bcFlagUni_->at(i) = 1;
                        }
                        //out
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0] + length - tol) &&
                            this->pointsUni_->at(i).at(1) > (coorRec[1] + tol) &&
                            this->pointsUni_->at(i).at(1) < (coorRec[1] + width - tol)&&
                            this->pointsUni_->at(i).at(2) > (coorRec[2] + tol) &&
                            this->pointsUni_->at(i).at(2) < (coorRec[2] + height - tol)) {
                            this->bcFlagUni_->at(i) = 3;
                        }
                    }
                    for (int i=0; i<this->pointsUni_->size(); i++) {
                        if (this->pointsRep_->at(i).at(0) < (coorRec[0] + tol) ) {
                            this->bcFlagRep_->at(i) = 2;
                        }
                        //bottom
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0] + tol) &&
                            this->pointsRep_->at(i).at(2) < (coorRec[2] + tol) ) {
                            this->bcFlagRep_->at(i) = 1;
                        }
                        //top
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0] + tol) &&
                            this->pointsRep_->at(i).at(2) > (coorRec[2] + height - tol) ) {
                            this->bcFlagRep_->at(i) = 1;
                        }
                        //front
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0] + tol) &&
                            this->pointsRep_->at(i).at(1) < (coorRec[1] + tol) ) {
                            this->bcFlagRep_->at(i) = 1;
                        }
                        //back
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0] + tol) &&
                            this->pointsRep_->at(i).at(1) > (coorRec[1] + width - tol) ) {
                            this->bcFlagRep_->at(i) = 1;
                        }
                        //out
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0] + length - tol) &&
                            this->pointsRep_->at(i).at(1) > (coorRec[1] + tol) &&
                            this->pointsRep_->at(i).at(1) < (coorRec[1] + width - tol)&&
                            this->pointsRep_->at(i).at(2) > (coorRec[2] + tol) &&
                            this->pointsRep_->at(i).at(2) < (coorRec[2] + height - tol)) {
                            this->bcFlagRep_->at(i) = 3;
                        }
                    }
                    break;
                case 2: //BFS
                    for (int i=0; i<this->pointsUni_->size(); i++) {
                        if (this->pointsUni_->at(i).at(0) < (coorRec[0] +tol) && this->pointsUni_->at(i).at(1) > (coorRec[1] - tol) && this->pointsUni_->at(i).at(1) < (coorRec[1]+ width +tol)
                            && this->pointsUni_->at(i).at(2) > (coorRec[2]+1. -tol) && this->pointsUni_->at(i).at(1) < (coorRec[2]+ height +tol)) {
                            this->bcFlagUni_->at(i) = 2;
                        }
                        //bottom top
                        if (this->pointsUni_->at(i).at(2) < (coorRec[2] + tol) || (this->pointsUni_->at(i).at(2) < (coorRec[2]+1. + tol) && this->pointsUni_->at(i).at(0) < (coorRec[0]+1. + tol)) || this->pointsUni_->at(i).at(2) > (coorRec[2]+height - tol)) {
                            this->bcFlagUni_->at(i) = 1;
                        }
                        //front back
                        if (this->pointsUni_->at(i).at(1) < (coorRec[1] + tol) || this->pointsUni_->at(i).at(1) > (coorRec[1]+width - tol)) {
                            this->bcFlagUni_->at(i) = 1;
                        }
                        //step left
                        if (this->pointsUni_->at(i).at(0) < (coorRec[0]+1. + tol) && this->pointsUni_->at(i).at(0) > (coorRec[0]+1. - tol) && this->pointsUni_->at(i).at(2) < (coorRec[2]+1.+tol)) {
                            this->bcFlagUni_->at(i) = 1;
                        }
                        if (this->pointsUni_->at(i).at(0) > (coorRec[0]+length - tol) && this->pointsUni_->at(i).at(1) > (coorRec[1] + tol) && this->pointsUni_->at(i).at(1) < (coorRec[1] + width - tol)
                            && this->pointsUni_->at(i).at(2) > (coorRec[2] + tol) && this->pointsUni_->at(i).at(2) < (coorRec[2] + height - tol)) {
                            this->bcFlagUni_->at(i) = 3;
                        }
                    }
                    for (int i=0; i<this->pointsRep_->size(); i++) {
                        if (this->pointsRep_->at(i).at(0) < (coorRec[0] +tol) && this->pointsRep_->at(i).at(1) > (coorRec[1] - tol) && this->pointsRep_->at(i).at(1) < (coorRec[1]+ width +tol)
                            && this->pointsRep_->at(i).at(2) > (coorRec[2]+1. -tol) && this->pointsRep_->at(i).at(1) < (coorRec[2]+ height +tol)) {
                            this->bcFlagRep_->at(i) = 2;
                        }
                        //bottom top
                        if (this->pointsRep_->at(i).at(2) < (coorRec[2] + tol) || (this->pointsRep_->at(i).at(2) < (coorRec[2]+1. + tol) && this->pointsRep_->at(i).at(0) < (coorRec[0]+1. + tol)) || this->pointsRep_->at(i).at(2) > (coorRec[2]+height - tol)) {
                            this->bcFlagRep_->at(i) = 1;
                        }
                        //front back
                        if (this->pointsRep_->at(i).at(1) < (coorRec[1] + tol) || this->pointsRep_->at(i).at(1) > (coorRec[1]+width - tol)) {
                            this->bcFlagRep_->at(i) = 1;
                        }
                        //step left
                        if (this->pointsRep_->at(i).at(0) < (coorRec[0]+1. + tol) && this->pointsRep_->at(i).at(0) > (coorRec[0]+1. - tol) && this->pointsRep_->at(i).at(2) < (coorRec[2]+1.+tol)) {
                            this->bcFlagRep_->at(i) = 1;
                        }
                        if (this->pointsRep_->at(i).at(0) > (coorRec[0]+length - tol) && this->pointsRep_->at(i).at(1) > (coorRec[1] + tol) && this->pointsRep_->at(i).at(1) < (coorRec[1] + width - tol)
                            && this->pointsRep_->at(i).at(2) > (coorRec[2] + tol) && this->pointsRep_->at(i).at(2) < (coorRec[2] + height - tol)) {
                            this->bcFlagRep_->at(i) = 3;
                        }
                    }
                    break;
                default:
                    break;
            }
            
            break;
        default:
            break;
    }
}
    
template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::buildElementMap(){
    
    Teuchos::Array<GO> elementsGlobalMapping( this->elementsC_->numberElements() );
    LO offset = this->comm_->getRank() * elementsGlobalMapping.size();
    for (int i=0; i<elementsGlobalMapping.size(); i++)
        elementsGlobalMapping[i] = i + offset;
    
    std::string underlyingLib = this->mapRepeated_->getUnderlyingLib();
    this->elementMap_.reset(new Map<LO,GO,NO>( underlyingLib, (GO) -1, elementsGlobalMapping(), 0, this->comm_) );
    
}

// We assemble the surface elements to the corresponding mesh
// In 2D we build edgeElements and in 3D the triangle Elements

/* template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::assembleSurfaceElements(){




}

template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::assignEdgeFlags(){

	vec2D_LO_Type markedPoints(0);
	elements = this->elements_;
	edgeElements = this->edgeElements_;
	MapConstPtr_Type mapRepeatedP1 = this->getMapRepeated();
	vec_int_Type newFlags(edgeElements->numberElements());

	MapConstPtr_Type edgeMap = this->getEdgeMap();

	vec_int_Type markedTrue(edgeElements->numberElements());
	for(int i=0; i< edgeElements->numberElements() ; i++){
		LO p1ID =edgeElements->getElement(i).getNode(0);
		LO p2ID =edgeElements->getElement(i).getNode(1);
		newFlags[i]=meshP1->determineFlagP2( meshP1, p1ID, p2ID, i , markedPoints );
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
			mapGlobalInterfaceUnique = mapGlobalInterface->buildUniqueMap( meshP1->rankRange_ );
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

template <class SC, class LO, class GO, class NO>
void MeshStructured<SC,LO,GO,NO>::buildSurfaceTriangleElements(ElementsPtr_Type elements, EdgeElementsPtr_Type edgeElements, SurfaceElementsPtr_Type surfaceTriangleElements, MapConstPtr_Type edgeMap, MapConstPtr_Type elementMap ){
    TEUCHOS_TEST_FOR_EXCEPTION( elements.is_null(), std::runtime_error, "Elements not initialized!");
    TEUCHOS_TEST_FOR_EXCEPTION( surfaceTriangleElements.is_null(), std::runtime_error, "Surface Triangle Elements not initialized!");


	vec_LO_Type nodeInd(4);
	vec2D_int_Type newTriangles(4,vec_int_Type(0)); // new Triangles
	
	vec_GO_Type globalInterfaceIDs = edgeElements->determineInterfaceEdges(edgeMap);
	//if(edgeElements->getEdgesOfElement(0) ) here we need some sort of test if the function was already used
	edgeElements->matchEdgesToElements(elementMap);


	for(int i=0; i<elements->numberElements(); i++){
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(i); // indeces of edges belonging to element

		// Extract the four points of tetraeder
		vec_int_Type nodeInd(0);
		for(int i=0; i<6; i++)	{
			nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(0));
			nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(1));
		}
		sort( nodeInd.begin(), nodeInd.end() );
		nodeInd.erase( unique( nodeInd.begin(), nodeInd.end() ), nodeInd.end() );

		// With our sorted Nodes we construct edges as follows

		// Edge_0 = [x_0,x_1]
		// Edge_1 = [x_0,x_2]	 
		// Edge_2 = [x_0,x_3]	 
		// Edge_3 = [x_1,x_2]
		// Edge_4 = [x_1,x_3]	 
		// Edge_5 = [x_2,x_3]

		vec_int_Type edgeNumbersTmp = edgeNumbers;
		for(int i=0; i<6; i++){
			if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[0] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[0]){
				if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[1] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[1])
					edgeNumbers[0] = edgeNumbersTmp[i];
				else if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[2] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[2])
					edgeNumbers[1] = edgeNumbersTmp[i];
				else if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[3] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[3])
					edgeNumbers[2] = edgeNumbersTmp[i];	
			}
			else if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[1] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[1]){
				if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[2] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[2])
					edgeNumbers[3] = edgeNumbersTmp[i];
				else if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[3] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[3])
					edgeNumbers[4] = edgeNumbersTmp[i];
			}
			 else if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[2] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[2]){
				if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[3] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[3])
					edgeNumbers[5] = edgeNumbersTmp[i];
			}
		}

		// We hace 4 Triangles in our Tetraedron
		// If one or more of those Triangles are Part of the domains' boundaries, they are added to the element in question as subelement
		// We extract them in follwing pattern:

		// Tri_0 = [x_0,x_1,x_2] -> Edge 0,1,3
		// Tri_1 = [x_0,x_1,x_3] -> Edge 0,2,4
		// Tri_2 = [x_0,x_2,x_3] -> Edge 1,2,5
		// Tri_3 = [x_1,x_2,x_3] -> Edge 3,4,5

		// We check if one or more of these triangles are part of the boundary surface and determine there flag

		vec2D_int_Type originTriangles(4,vec_int_Type(3));
		originTriangles[0] = {nodeInd[0],nodeInd[1],nodeInd[2]};
		originTriangles[1] = {nodeInd[0],nodeInd[1],nodeInd[3]};
		originTriangles[2] = {nodeInd[0],nodeInd[2],nodeInd[3]};
		originTriangles[3] = {nodeInd[1],nodeInd[2],nodeInd[3]};
		
		
		vec_int_Type originFlag(4,10); // Triangle Flag

		int numberSubElSurf=0;
		vec_LO_Type triTmp(3);
		vec_int_Type originTriangleTmp(3);
		int entry; 
		if (elements->getElement(i).subElementsInitialized() ){
			numberSubElSurf = elements->getElement(i).getSubElements()->numberElements();
			for(int k=0; k< numberSubElSurf ; k++){
				triTmp =elements->getElement(i).getSubElements()->getElement(k).getVectorNodeList();
				for(int j=0; j<4 ; j++){
					originTriangleTmp = originTriangles[j];
					sort(originTriangleTmp.begin(),originTriangleTmp.end());
					sort(triTmp.begin(),triTmp.end());
					if(triTmp[0] == originTriangleTmp[0] && triTmp[1] == originTriangleTmp[1] &&  triTmp[2] == originTriangleTmp[2] ) 
						originFlag[j] = elements->getElement(i).getSubElements()->getElement(k).getFlag();
				
				}
			}
		}

		// Furthermore we have to determine whether the triangles are part of the interface between processors, as we need this information to determine if edges
		// that emerge on the triangles are part of the interface
		// A triangle is part of the interface if all of its edges are part of the interface (the information if edges are part of the interface was determined
		// in the beginning of the Mesh Refinement by 'determineInterfaceEdges')

		vec_bool_Type interfaceSurface = checkInterfaceSurface(edgeElements,originFlag, edgeNumbers,i);
		
		for(int j=0; j<4; j++){	
			sort( newTriangles.at(j).begin(), newTriangles.at(j).end() );
			FiniteElement feNew(originTriangles[j],originFlag[j]);
			feNew.setInterfaceElement(interfaceSurface[j]);
			surfaceTriangleElements->addSurface(feNew, i);
		}
	}
	vec2D_GO_Type combinedSurfaceElements;
	surfaceTriangleElements->sortUniqueAndSetGlobalIDsParallel(elementMap,combinedSurfaceElements);
	
	surfaceTriangleElements->setElementsSurface( combinedSurfaceElements );

	surfaceTriangleElements->setUpElementsOfSurface(elementMap,edgeMap, edgeElements);	
}
  */  
}
#endif
