#ifndef AssembleFE_NeoHooke2D_DEF_hpp
#define AssembleFE_NeoHooke2D_DEF_hpp

#include "AssembleFE_NeoHooke2D_decl.hpp"
#ifdef FEDD_HAVE_ACEGENINTERFACE
#include "aceinterface.hpp"
#endif

#include <vector>
// #include <iostream>

namespace FEDD
{

	template <class SC, class LO, class GO, class NO>
	AssembleFE_NeoHooke2D<SC, LO, GO, NO>::AssembleFE_NeoHooke2D(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params, tuple_disk_vec_ptr_Type tuple) : AssembleFE<SC, LO, GO, NO>(flag, nodesRefConfig, params, tuple)
	{
		// Extracting values from ParameterList
		int numMaterials =  this->params_->sublist("Parameter Solid").get("Number of Materials", 1);
		int materialID = 1;

		E_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("E", 0.38);
		poissonRatio_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("Poisson Ratio", 0.49e-0);
		bodyForceX_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("Body Force X", 0.1e-0);
		bodyForceY_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("Body Force Y", 0.0e-0);
		rho_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("Density",1.0);
		t_ =  this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("Thickness",1.0);

		//cout << "--- Init AssembleFE_NeoHooke2D Element --- EMOD " << E_  << " Body Force X " << bodyForceX_ << endl;

		FEType_ = std::get<1>(this->diskTuple_->at(0));	   // FEType of Disk
		dofs_ = std::get<2>(this->diskTuple_->at(0)); // Degrees of freedom per node
		numNodes_ = std::get<3>(this->diskTuple_->at(0)); // Number of nodes of element

		dofsElement_ = dofs_ * numNodes_ ; // "Dimension of return matrix"

	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_NeoHooke2D<SC, LO, GO, NO>::assembleJacobian()
	{

		SmallMatrixPtr_Type elementMatrix = Teuchos::rcp(new SmallMatrix_Type(dofsElement_, 0.));

		assembleNeoHook2D(elementMatrix);
		// elementMatrix->print();
		this->jacobian_ = elementMatrix;
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_NeoHooke2D<SC, LO, GO, NO>::assembleRHS()
	{

		this->rhsVec_.reset(new vec_dbl_Type(dofsElement_, 0.));
#ifdef FEDD_HAVE_ACEGENINTERFACE

		std::vector<double> positions(dofsElement_);
		int count = 0;
		for (int i = 0; i < numNodes_; i++)
		{
			for (int j = 0; j < dofs_; j++)
			{
				positions[count] = this->getNodesRefConfig()[i][j];
				count++;
			}
		}

		std::vector<double> displacements(dofsElement_);
		for (int i = 0; i < dofsElement_; i++){
			displacements[i] = (*this->solution_)[i];
		}

		std::vector<double> domainData(6);
		domainData[0] = E_;
		domainData[1] = poissonRatio_;
		domainData[2] = bodyForceX_;
		domainData[3] = bodyForceY_;
		domainData[4] = rho_;
		domainData[5] = t_;

		AceGenInterface::NeoHookeTriangle2D3PlaneStress elem(&positions[0], &displacements[0], domainData, this->globalElementID_);
		elem.computeResiduum();
		
		std::vector<double> residuum = elem.getResiduum();

		for (int i = 0; i < dofsElement_; i++){
			(*this->rhsVec_)[i] = -residuum[i];
			cout << (*this->rhsVec_)[i] << " " ;
		}
		cout << endl;


#endif
	
	}


	template <class SC, class LO, class GO, class NO>
	void AssembleFE_NeoHooke2D<SC, LO, GO, NO>::assembleNeoHook2D(SmallMatrixPtr_Type &elementMatrix)
	{
#ifdef FEDD_HAVE_ACEGENINTERFACE

		std::vector<double> positions(dofsElement_);
		int count = 0;
		cout << "--------- Positions ----------" << endl;
		for (int i = 0; i < numNodes_; i++)
		{
			for (int j = 0; j < dofs_; j++)
			{
				positions[count] = this->getNodesRefConfig()[i][j];
				cout << positions[count] << " " ;
				count++;
			}
		}
		cout << endl;

		std::vector<double> displacements(dofsElement_);
		cout << "--------- Displacements ----------" << endl;
		for (int i = 0; i < dofsElement_; i++){
			displacements[i] = (*this->solution_)[i];
			cout << displacements[i] << " " ;
		}
		cout << endl;

		std::vector<double> domainData(6);
		domainData[0] = E_;
		domainData[1] = poissonRatio_;
		domainData[2] = bodyForceX_;
		domainData[3] = bodyForceY_;
		domainData[4] = rho_;
		domainData[5] = t_;
		cout << " Domain Data " << domainData[0] << " " << domainData[1] << " " <<domainData[2] << " " <<domainData[3] << " " <<domainData[4] << " " <<domainData[5]  << endl;
		AceGenInterface::NeoHookeTriangle2D3PlaneStress elem = AceGenInterface::NeoHookeTriangle2D3PlaneStress(&positions[0], &displacements[0], domainData, this->globalElementID_);
		elem.compute();

		double **stiffnessMatrix = elem.getStiffnessMatrix();
		

		for (UN i = 0; i < this->dofsElement_; i++)
		{
			for (UN j = 0; j < this->dofsElement_; j++)
			{
				(*elementMatrix)[i][j] = stiffnessMatrix[i][j];
			}
		}
		if(globalElementID_==0)
			elementMatrix->print();
	
	#endif
	}
} // namespace FEDD
#endif // AssembleFE_NeoHooke2D_DEF_hpp
