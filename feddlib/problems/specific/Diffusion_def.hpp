#ifndef DIFFUSION_def_hpp
#define DIFFUSION_def_hpp
#include "Diffusion_decl.hpp"
/*!
 Definition of Diffusion
 
 @brief Diffusion
 @author Lea Sassmannshausen
 @version 1.0
 @copyright LS
 */

namespace FEDD {

template<class SC,class LO,class GO,class NO>
Diffusion<SC,LO,GO,NO>::Diffusion(const DomainConstPtr_Type &domain, std::string FEType, ParameterListPtr_Type parameterList, vec2D_dbl_Type diffusionTensor, bool vectorDiffusion):
Problem<SC,LO,GO,NO>(parameterList, domain->getComm())
//vectorDiffusion_(vectorDiffusion)
{
 
    this->addVariable( domain , FEType , "u" , 1);
    this->dim_ = this->getDomain(0)->getDimension();
	this->diffusionTensor_ = diffusionTensor;
}

/*template<class SC,class LO,class GO,class NO>
Diffusion<SC,LO,GO,NO>::~Diffusion(){

}*/
    
template<class SC,class LO,class GO,class NO>
void Diffusion<SC,LO,GO,NO>::info(){
    this->infoProblem();
}

template<class SC,class LO,class GO,class NO>
void Diffusion<SC,LO,GO,NO>::assemble( std::string type ) const{
    
    if (this->verbose_)
        std::cout << "-- Assembly Laplace with Diffusion Tensor ... " << std::flush;

    MatrixPtr_Type A;
    vec_dbl_Type funcParameter(1,0.);


    A = Teuchos::rcp(new Matrix_Type( this->domainPtr_vec_.at(0)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
    this->feFactory_->assemblyLaplaceDiffusion(this->dim_, this->domain_FEType_vec_.at(0), 2, A, this->diffusionTensor_ );


    this->system_->addBlock(A,0,0);
    
    this->assembleSourceTerm( 0. );
    this->addToRhs( this->sourceTerm_ );
    
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

template<class SC,class LO,class GO,class NO>
typename Diffusion<SC,LO,GO,NO>::MatrixPtr_Type Diffusion<SC,LO,GO,NO>::getMassMatrix() const{
	
    MatrixPtr_Type A;
	A = Teuchos::rcp(new Matrix_Type( this->domainPtr_vec_.at(0)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
	this->feFactory_->assemblyMass(this->dim_,this->domain_FEType_vec_.at(0),"Scalar", A);

	return A;

}
    
}
#endif
