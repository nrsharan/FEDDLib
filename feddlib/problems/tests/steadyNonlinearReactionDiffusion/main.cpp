#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/DiffusionReaction.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include "feddlib/problems/Solver/DAESolverInTime.hpp"


/*!
 main of unsteadyLinearDiffusion problem

 @brief unsteadyLinearDiffusion main
 @author Lea Saßmannshausen
 @version 1.0
 @copyright LS
 */

void zeroBC(double* x, double* res, double t, const double* parameters){
    res[0] = 0.;
}
void oneBC(double* x, double* res, double t, const double* parameters){
    res[0] = 1.;
}
void twoBC(double* x, double* res, double t, const double* parameters){
    res[0] = 2.;
}
void threeBC(double* x, double* res, double t, const double* parameters){
    res[0] = 3.;
}
void zeroBC2D(double* x, double* res, double t, const double* parameters){
    res[0] = 0.;
    res[1] = 0.;
}
void zeroBC3D(double* x, double* res, double t, const double* parameters){
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
}

void oneFunc(double* x, double* res, double* parameters){
    res[0] = 1.;
}

void zeroFunc(double* x, double* res, double* parameters){
    res[0] = 1.;
}

void reactionFunc(double* x, double* res, double* parameters){
	double E0 = parameters[0];
	double E1 = parameters[1];
    res[0] = E0 - (E0-E1)* ( x[0] / (x[0] + 1./2) ) ;

}

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;

int main(int argc, char *argv[]) {
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
	typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;

    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");


    std::string vectorLaplace = "false";
    myCLP.setOption("vectorLaplace",&vectorLaplace,"vectorLaplace");
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");
    double length = 4.;
    myCLP.setOption("length",&length,"length of domain.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    bool vL ( !vectorLaplace.compare("true") );
    
    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        parameterListAll->setParameters(*parameterListPrec);
        parameterListAll->setParameters(*parameterListSolver);

        // Mesh
        int dim = parameterListProblem->sublist("Parameter").get("Dimension",2);
        int m = parameterListProblem->sublist("Parameter").get("H/h",5);
        std::string FEType = parameterListProblem->sublist("Parameter").get("Discretization","P1");
        std::string meshType = parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        std::string meshName = parameterListProblem->sublist("Parameter").get("Mesh Name","");
        std::string meshDelimiter = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");

        int n;
        int size = comm->getSize();
        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        size -= numProcsCoarseSolve;

        int minNumberSubdomains;
        if (!meshType.compare("structured") || !meshType.compare("unstructured_struct")) {
            minNumberSubdomains = 1;
        }
        else if(!meshType.compare("structured_bfs") || !meshType.compare("unstructured_bfs")){
            minNumberSubdomains = (int) 2*length+1;
        }


        Teuchos::RCP<Domain<SC,LO,GO,NO> > domain;

        Teuchos::RCP<Domain<SC,LO,GO,NO> > domainP1;
        Teuchos::RCP<Domain<SC,LO,GO,NO> > domainP2;
        domainP1.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
        
        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
        domainP1Array[0] = domainP1;
        
        ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );
        MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
        
        partitionerP1.readAndPartition();

        if (FEType=="P2") {
            domainP2.reset( new Domain<SC,LO,GO,NO>( comm, dim ));
            domainP2->buildP2ofP1Domain( domainP1 );
            domain = domainP2;
        }
        else
            domain = domainP1;
        



        // ####################
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory(new BCBuilder<SC,LO,GO,NO>( ));
        if (vL) {
            if (dim==2) {
                bcFactory->addBC(zeroBC2D, 1, 0, domain, "Dirichlet", dim);
                bcFactory->addBC(zeroBC2D, 2, 0, domain, "Dirichlet", dim);
                bcFactory->addBC(zeroBC2D, 3, 0, domain, "Dirichlet", dim);
            }
            else if(dim==3){
                bcFactory->addBC(zeroBC3D, 1, 0, domain, "Dirichlet", dim);
                bcFactory->addBC(zeroBC3D, 2, 0, domain, "Dirichlet", dim);
                bcFactory->addBC(zeroBC3D, 3, 0, domain, "Dirichlet", dim);
            }
        }
        else{
            bcFactory->addBC(threeBC, 1, 0, domain, "Dirichlet", 1);
            bcFactory->addBC(zeroBC, 2, 0, domain, "Dirichlet", 1);
            bcFactory->addBC(zeroBC, 3, 0, domain, "Dirichlet", 1);
            bcFactory->addBC(zeroBC, 4, 0, domain, "Dirichlet", 1);
        }

		MultiVectorPtr_Type diffusionTensor = Teuchos::RCP<MultiVector_Type>(new MultiVector_Type ( domain->getElementMap(),1));
		diffusionTensor->putScalar(1.);
		Teuchos::ArrayRCP< SC >  linearDiff = diffusionTensor->getDataNonConst( 0 );
		for(int j=0; j<domain->getElementsC()->numberElements(); j++){
			for(int k = 0; k< 3 ; k++){
				int i = domain->getElementsC()->getElement(j).getNode(k);
				if(domain->getPointsRepeated()->at(i).at(0) == 0.5 &&  domain->getPointsRepeated()->at(i).at(1)>0 && domain->getPointsRepeated()->at(i).at(1) <1)
					linearDiff[j] = 100000.;
			}
		}
        DiffusionReaction<SC,LO,GO,NO> diffusionReaction(domain,FEType,parameterListAll, diffusionTensor,reactionFunc);
        {
       
            domain->info();
            diffusionReaction.info();
            
            {

                diffusionReaction.addBoundaries( bcFactory );

                diffusionReaction.addRhsFunction( zeroFunc );
     
                diffusionReaction.initializeProblem();
                diffusionReaction.assemble();
                                    
                diffusionReaction.setBoundaries();
                
                std::string nlSolverType = parameterListProblem->sublist("General").get("Linearization","Newton");
                NonLinearSolver<SC,LO,GO,NO> nlSolver(nlSolverType);
                nlSolver.solve( diffusionReaction );

            }


       		Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolution = diffusionReaction.getSolution()->getBlock(0);

            exPara->setup("solutionDiffusionReaction", domain->getMesh(), FEType);
            
            exPara->addVariable(exportSolution, "u", "Scalar", 1, domain->getMapUnique(), domain->getMapUniqueP2());

            exPara->save(0.0);

         }





        

       
    }
    return(EXIT_SUCCESS);
}
