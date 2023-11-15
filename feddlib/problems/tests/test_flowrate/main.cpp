#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/FSI.hpp"
#include "feddlib/problems/specific/Laplace.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*! Test case for specific artery geometrie or straight tube geometry. Inflow depends on inflow region
	-> artery: Inflow scaled with normal vector on inflow (x,y,z) * laplaceInflow	

*/



void zeroBC(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void parabolicInflow3D(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired flow rate
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing

    res[0] = parameters[3];

    return;
}

void parabolicInflowDirection3D(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing

    res[0] = 0.;
    res[1] = 0.;
    res[2] = parameters[0] * x[0];

    return;
}


void rhsDummy(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void oneFunc(double* x, double* res, double* parameters){
    res[0] = 1.;
}

void dummyFunc(double* x, double* res, double t, const double* parameters)
{
    return;
}

void reactionFunc(double* x, double* res, double* parameters){
	
    double m = 0.0;	
    res[0] = m * x[0];

}


typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;

using namespace FEDD;
using namespace Teuchos;
using namespace std;

int main(int argc, char *argv[])
{


    typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef Mesh<SC,LO,GO,NO> Mesh_Type;
    typedef RCP<Mesh_Type> MeshPtr_Type;
    typedef RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef ExporterParaView<SC,LO,GO,NO> ExporterPV_Type;
    typedef RCP<ExporterPV_Type> ExporterPVPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

  
    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string xmlProbL = "plistProblemLaplace.xml";
    myCLP.setOption("probLaplace",&xmlProbL,".xml file with Inputparameters.");
    string xmlPrecL = "plistPrecLaplace.xml";
    myCLP.setOption("precLaplace",&xmlPrecL,".xml file with Inputparameters.");
    string xmlSolverL = "plistSolverLaplace.xml";
    myCLP.setOption("solverLaplace",&xmlSolverL,".xml file with Inputparameters.");
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    bool verbose (comm->getRank() == 0);

    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);   
        
        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",2);
        
        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","unstructured");
        
        string      discType        = parameterListProblem->sublist("Parameter").get("Discretization","P2");
        string preconditionerMethod = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        int         n;

        TimePtr_Type totalTime(TimeMonitor_Type::getNewCounter("FEDD - main - Total Time"));
        TimePtr_Type buildMesh(TimeMonitor_Type::getNewCounter("FEDD - main - Build Mesh"));

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);

        int size = comm->getSize() - numProcsCoarseSolve;

        // #####################
        // Mesh bauen und wahlen
        // #####################
        if (verbose)
        {
            cout << "###############################################" <<endl;
            cout << "############ Starting test  ... ################" <<endl;
            cout << "###############################################" <<endl;
        }

        DomainPtr_Type domainP1fluid;      
       	DomainPtr_Type domainP2fluid;
       
        DomainPtr_Type domainFluidVelocity;
        DomainPtr_Type domainFluidPressure;
        
        TimeMonitor_Type totalTimeMonitor(*totalTime);
    
        TimeMonitor_Type buildMeshMonitor(*buildMesh);
        if (verbose)
        {
            cout << " -- Building Mesh ... " << flush;
        }

        domainP1fluid.reset( new Domain_Type( comm, dim ) );
   
        domainP2fluid.reset( new Domain_Type( comm, dim ) );
                            
                        
        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
        domainP1Array[0] = domainP1fluid;

        
        ParameterListPtr_Type pListPartitioner = sublist( parameterListProblem, "Mesh Partitioner" );

        pListPartitioner->set("Build Edge List",true);
        pListPartitioner->set("Build Surface List",true);
        
        MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
        
        partitionerP1.readAndPartition(15);
        
        if (!discType.compare("P2")){
            domainP2fluid->buildP2ofP1Domain( domainP1fluid );
		}

        if (!discType.compare("P2"))
        {
            domainFluidVelocity = domainP2fluid;
            domainFluidPressure = domainP1fluid;
        }
        else
        {
            domainFluidVelocity = domainP1fluid;
            domainFluidPressure = domainP1fluid;

        }

        if (verbose){
            cout << "done! -- " << endl;
        }
            
                        

        // #####################
        // Problem definieren
        // #####################

            
        std::vector<double> parameter_vec(1, parameterListProblem->sublist("Parameter").get("Max Velocity",1.));
        parameter_vec.push_back( parameterListProblem->sublist("Parameter").get("Max Ramp Time",2.) );
            
        //#############################################
        //#############################################
        //#### Compute parabolic inflow with laplacian
        //#############################################
        //#############################################
        MultiVectorPtr_Type solutionLaplace;
        MultiVectorConstPtr_Type solutionLaplaceConst;
        {
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryLaplace(new BCBuilder<SC,LO,GO,NO>( ));
            
            bcFactoryLaplace->addBC(zeroBC, 2, 0, domainFluidVelocity, "Dirichlet", 1); //inflow ring
            bcFactoryLaplace->addBC(zeroBC, 3, 0, domainFluidVelocity, "Dirichlet", 1); //outflow ring
            bcFactoryLaplace->addBC(zeroBC, 6, 0, domainFluidVelocity, "Dirichlet", 1); //surface
            
            ParameterListPtr_Type parameterListProblemL = Teuchos::getParametersFromXmlFile(xmlProbL);
            ParameterListPtr_Type parameterListPrecL = Teuchos::getParametersFromXmlFile(xmlPrecL);
            ParameterListPtr_Type parameterListSolverL = Teuchos::getParametersFromXmlFile(xmlSolverL);

            ParameterListPtr_Type parameterListLaplace(new Teuchos::ParameterList(*parameterListProblemL)) ;
            parameterListLaplace->setParameters(*parameterListPrecL);
            parameterListLaplace->setParameters(*parameterListSolverL);
            
            Laplace<SC,LO,GO,NO> laplace( domainFluidVelocity, discType, parameterListLaplace, false );
            {
                laplace.addRhsFunction(oneFunc);
                laplace.addBoundaries(bcFactoryLaplace);
                
                laplace.initializeProblem();
                laplace.assemble();
                laplace.setBoundaries();
                laplace.solve();
            }
            
            //We need the values in the inflow area. Therefore, we use the above bcFactory and the volume flag 10 and the outlet flag 5 and set zero Dirichlet boundary values
            bcFactoryLaplace->addBC(zeroBC, 5, 0, domainFluidVelocity, "Dirichlet", 1);
            bcFactoryLaplace->addBC(zeroBC, 15, 0, domainFluidVelocity, "Dirichlet", 1);
            bcFactoryLaplace->setRHS( laplace.getSolution(), 0.);
            solutionLaplace = laplace.getSolution()->getBlockNonConst(0);
        
            SC maxValue = solutionLaplace->getMax();
            solutionLaplace->scale(1./maxValue); // normalizing solution
            
            solutionLaplaceConst = solutionLaplace; 

            parameter_vec.push_back(1.0); // We scaled the solution beforehand, so we dont need the actual maxValue any more and replace it with 1.


            parameter_vec.push_back(parameterListProblem->sublist("Parameter").get("Flowrate",3.0));

            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
            
            exPara->setup("parabolicInflow", domainFluidVelocity->getMesh(), discType);
  
            MultiVectorConstPtr_Type valuesConst = laplace.getSolution()->getBlock(0);
            exPara->addVariable( valuesConst, "values", "Scalar", 1, domainFluidVelocity->getMapUnique() );

            exPara->save(0.0);
            exPara->closeExporter();

        }
        
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );


        // Fluid-RW
        {

            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryFluid( new BCBuilder<SC,LO,GO,NO>( ) );
                            
            bcFactory->addBC(parabolicInflow3D, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplaceConst, true , parabolicInflowDirection3D); // inflow 


            bcFactory->addBC(parabolicInflow3D, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplaceConst, true , parabolicInflowDirection3D); // inflow 
 
        }
        std::vector<MapConstPtr_Type> maps(1);
        maps[0] = domainFluidVelocity->getMapVecFieldUnique();
		BlockMultiVectorPtr_Type  vector = Teuchos::rcp( new BlockMultiVector_Type( maps ) ) ;
        bcFactory->setRHS(vector);
    	
		
		FE<SC,LO,GO,NO> fe;
		fe.addFE(domainFluidVelocity);
		
		double flowRateInlet;
		double flowRate = parameterListProblem->sublist("Parameter").get("Flowrate",3.0);
		MultiVectorPtr_Type vector_rep = Teuchos::rcp(new MultiVector_Type ( domainFluidVelocity->getMapVecFieldRepeated() ) );   
    	vector_rep->importFromVector(vector->getBlock(0),false,"Insert");
    	
		fe.assemblyFlowRate(dim, flowRateInlet, discType , dim, 4 , vector_rep);
		
		
		cout << " Flowrate set by parameterlist: " << flowRate << "  |  flowrate set at inlet: " << flowRateInlet << endl;
    	TEUCHOS_TEST_FOR_EXCEPTION( abs(flowRateInlet - flowRate) > 1e-12 , std::runtime_error, " Flowrate meassured != flowrate prescribed");

    
    }
    

    TimeMonitor_Type::report(std::cout);

    return(EXIT_SUCCESS);
}
