#ifndef MAIN_TIMER_START
#define MAIN_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Main") + std::string(S))));
#endif

#ifndef MAIN_TIMER_STOP
#define MAIN_TIMER_STOP(A) A.reset();
#endif

#include "feddlib/amr/AdaptiveMeshRefinement.hpp"

#include <Teuchos_TestForException.hpp>

#include "feddlib/problems/abstract/Problem.hpp"

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"

#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"

#include "feddlib/problems/Solver/DAESolverInTime.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 main of Stokes problem
 
 @brief Unsteady Navier-Stokes main
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


void rhs0( double* p, double* res, const double* parameters){

	res[0] =0;
    res[1] =0;
	res[2] =0;

	//cout << " res[0] " << res[0] << " res[1] " << res[1] << endl;
}

void one(double* x, double* res, double t, const double* parameters){
    
    res[0] = 1.;
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}
void two(double* x, double* res, double t, const double* parameters){
    
    res[0] = 2.;
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}
void three(double* x, double* res, double t, const double* parameters){
    
    res[0] = 3.;
    res[1] = 0.;
	res[2] = 0.;
    
    return;
}
void four(double* x, double* res, double t, const double* parameters){
    
    res[0] =parameters[0]*1.;
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}
void five(double* x, double* res, double t, const double* parameters){
    
    res[0] = parameters[0]*1.;
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}
void six(double* x, double* res, double t, const double* parameters){
    
    res[0] = 6.;
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}

void zeroDirichlet(double* x, double* res, double t, const double* parameters){
    
    res[0] = 0.;
    
    return;
}

void zeroDirichlet2D(double* x, double* res, double t, const double* parameters){
    
    res[0] = 0.;
    res[1] = 0.;
    
    return;
}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters){
    
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}

void inflowParabolic2D(double* x, double* res, double t, const double* parameters){
    
    double H = parameters[1];
    res[0] = 4*parameters[0]*x[1]*(H-x[1])/(H*H);
    res[1] = 0.;
    
    return;
}

void inflowParabolic3D(double* x, double* res, double t, const double* parameters){
    
    double H = parameters[1];
    res[0] = 16*parameters[0]*x[1]*(H-x[1])*x[2]*(H-x[2])/(H*H*H*H);
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}

void dummyFunc(double* x, double* res, double t, const double* parameters){
    
    return ;
}

void dummyFuncSol(double* x, double* res){
    
	res[0] = 0.;

    return;
}

void inflowParabolic2DSin(double* x, double* res, double t, const double* parameters){

    double H = parameters[1];
    res[0] = sin(M_PI*t*0.125)*( 6*x[1]*(H-x[1]) ) / (H*H);
    res[1] = 0.;

    return;
}

void inflow3DRichter(double* x, double* res, double t, const double* parameters)
{
    double H = parameters[1];
    
    if(t < 1.)
    {
        res[0] = 9./8 * parameters[0] *x[1]*(H-x[1])*(H*H-x[2]*x[2])/( H*H*(H/2.)*(H/2.) ) * ((1 - cos(2.0*M_PI*t))/2.0);
        res[1] = 0.;
        res[2] = 0.;
    }
    else
    {
        res[0] = 9./8 * parameters[0] *x[1]*(H-x[1])*(H*H-x[2]*x[2])/( H*H*(H/2.)*(H/2.) );
        res[1] = 0.;
        res[2] = 0.;
    }
    
    return;
}


typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;

int main(int argc, char *argv[]) {
    
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    typedef Teuchos::RCP<Domain<SC,LO,GO,NO> > DomainPtr_Type;

	typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef Teuchos::RCP<Problem_Type> ProblemPtr_Type;

	typedef boost::function<void(double* x, double* res, double t, const double* parameters)>   BCFunc_Type;  

    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
    
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    bool verbose (comm->getRank() == 0);

    if (verbose) {
        cout << "###############################################################" <<endl;
        cout << "##################### Unsteady Navier-Stokes ####################" <<endl;
        cout << "###############################################################" <<endl;
    }

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    string xmlTekoPrecFile = "parametersTeko.xml";
    myCLP.setOption("tekoprecfile",&xmlTekoPrecFile,".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");

    double length = 4.;
    myCLP.setOption("length",&length,"length of domain.");
    
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        MPI_Finalize();
        return 0;
    }
    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);

        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);
        
        ParameterListPtr_Type parameterListPrecTeko = Teuchos::getParametersFromXmlFile(xmlTekoPrecFile);
        
        
        int dim = parameterListProblem->sublist("Parameter").get("Dimension",3);
        
        std::string discVelocity = parameterListProblem->sublist("Parameter").get("Discretization Velocity","P2");
        std::string discPressure = parameterListProblem->sublist("Parameter").get("Discretization Pressure","P1");

        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        string		meshName    	= parameterListProblem->sublist("Parameter").get("Mesh Name","some_mesh_here");
        string		meshDelimiter   = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");
        int 		volumeID        = parameterListProblem->sublist("Parameter").get("Volume ID",0);
        int 		m				= parameterListProblem->sublist("Parameter").get("H/h",5);        
        int         n;
        int			inflowID        = parameterListProblem->sublist("Parameter").get("Inflow ID",2);
        int			inflowBCID      = parameterListProblem->sublist("Parameter").get("InflowBC ID",-1);
        string      bcType          = parameterListProblem->sublist("Parameter").get("BC Type","parabolic");
        string      precMethod      = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
		int 		maxIter 		= parameterListProblem->sublist("Mesh Refinement").get("MaxIter",5);
		double maxVel 				= parameterListProblem->sublist("Parameter").get("MaxVelocity",2.);

		string modellProblem = parameterListProblem->sublist("Mesh Refinement").get("Modell Problem","Seminar1");

		std::vector<double> parameter_vec(0);
		parameter_vec.push_back(maxVel);//height of inflow region


		// Parameters determined by Modell Problem
		RhsFunc_Type rhs;
		Func_Type exactSolU;
		Func_Type exactSolP;
		BCFunc_Type flag1Func;
		BCFunc_Type flag2Func;
		BCFunc_Type flag3Func;
		BCFunc_Type flag4Func;
		BCFunc_Type flag5Func;
		BCFunc_Type flag6Func;
		BCFunc_Type flag7Func;


		if(modellProblem == "BFS" && dim ==2){
			rhs = rhs0;
			exactSolU = dummyFuncSol;
			exactSolP = dummyFuncSol;
			flag1Func = zeroDirichlet2D;
			flag2Func = inflowParabolic2D;
			parameter_vec.push_back(1);//height of inflow region	
		}
		else if(modellProblem == "BFS" && dim == 3){
			rhs = rhs0;
			exactSolU = dummyFuncSol;
			exactSolP = dummyFuncSol;
			flag1Func = zeroDirichlet3D;
			flag2Func = inflowParabolic3D;
			parameter_vec.push_back(1);//height of inflow region				
		}
		if(modellProblem == "Turek" && dim ==2){
			rhs = rhs0;
			exactSolU = dummyFuncSol;
			exactSolP = dummyFuncSol;
			flag1Func = zeroDirichlet2D;
			flag2Func = inflowParabolic2D;
			flag4Func = zeroDirichlet3D;
			parameter_vec.push_back(0.41);//height of inflow region	
		}
		else if(modellProblem == "Turek" && dim == 3){
			rhs = rhs0;
			exactSolU = dummyFuncSol;
			exactSolP = dummyFuncSol;
			flag1Func = zeroDirichlet3D;
			flag2Func = inflowParabolic3D;
			flag4Func = zeroDirichlet3D;
			parameter_vec.push_back(0.41);//height of inflow region				
		}
		
        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem));
        if (precMethod == "Monolithic")
            parameterListAll->setParameters(*parameterListPrec);
        else if(precMethod == "Teko")
            parameterListAll->setParameters(*parameterListPrecTeko);
        
        parameterListAll->setParameters(*parameterListSolver);
        
        int minNumberSubdomains;
        if (!meshType.compare("structured") || !meshType.compare("unstructured_struct")) {
            minNumberSubdomains = 1;
        }
        else if(!meshType.compare("structured_bfs") || !meshType.compare("unstructured_bfs")){
            minNumberSubdomains = (int) 2*length+1;
        }
        
        int size = comm->getSize();
        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        size -= numProcsCoarseSolve;
        int numProcsProblem = size;

        Teuchos::RCP<Teuchos::Time> totalTime(Teuchos::TimeMonitor::getNewCounter("main: Total Time"));
        Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
        Teuchos::RCP<Teuchos::Time> solveTime(Teuchos::TimeMonitor::getNewCounter("main: Solve problem time"));
        DomainPtr_Type domainPressure;
        DomainPtr_Type domainVelocity;
  
		domainPressure.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
		domainVelocity.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
		
		MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
		domainP1Array[0] = domainPressure;
		
		ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );
		MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
		
		partitionerP1.readAndPartition();
		
		Teuchos::RCP<Domain<SC,LO,GO,NO> > domainRefined;

		AdaptiveMeshRefinement<SC,LO,GO,NO> meshRefiner("NavierStokes",parameterListProblem,exactSolU,exactSolP); 
		
		

		int j=0;
		MAIN_TIMER_START(Total," Step 4:	 Total RefinementAlgorithm");
		//while(j<maxIter+1 ){

			MAIN_TIMER_START(buildP2," Step 0:	 buildP2Mesh");
			if (discVelocity=="P2" ) {
				domainVelocity.reset( new Domain<SC,LO,GO,NO>( comm, dim ));
		        domainVelocity->buildP2ofP1Domain( domainPressure );
				}
		    else
		        domainVelocity = domainPressure;
			
			MAIN_TIMER_STOP(buildP2);		

			MAIN_TIMER_START(Bounds," Step 1:	 bcFactory");
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );

			bcFactory->addBC(flag1Func, 1, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
			bcFactory->addBC(flag2Func, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
			//bcFactory->addBC(flag3Func, 3, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
			bcFactory->addBC(flag4Func, 4, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
			bcFactory->addBC(flag5Func, 5, 0, domainPressure, "Dirichlet", dim, parameter_vec);

      
			MAIN_TIMER_STOP(Bounds);	
			MAIN_TIMER_START(Solver," Step 2:	 solving PDE");

			//Teuchos::RCP<NavierStokes<SC,LO,GO,NO>> navierStokes( new NavierStokes<SC,LO,GO,NO> (domainVelocity, discVelocity, domainPressure, discPressure, parameterListAll ));

	        //navierStokes->info();

	        {
				//----------------------
				// Creating Navier Stokes 'specific problem type' object, adding boundary conditions, initializing problem, assembling and addin rhs boundaries
            	NavierStokes<SC,LO,GO,NO> navierStokes(domainVelocity, discVelocity, domainPressure, discPressure, parameterListAll );

	 			navierStokes.addBoundaries(bcFactory);
		        
		        navierStokes.initializeProblem();
		        
		        navierStokes.assemble();

		        navierStokes.setBoundariesRHS();

				// Differential Algebraic Equations solver. This contains the solving of different time steps and the corresponding solving of (non)linear systems
		        DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);
		        SmallMatrix<int> defTS(2);
		        defTS[0][0] = 1;
		        defTS[0][1] = 1;
		        defTS[1][0] = 0;
		        defTS[1][1] = 0;
				// Time stepping 
		        daeTimeSolver.defineTimeStepping(defTS);
				// Adding problem
		        daeTimeSolver.setProblem(navierStokes);

		        daeTimeSolver.setupTimeStepping();
				// Start of the time stepping process. We distinguish between (non)linear and Single-/Multi-step methods.
		        daeTimeSolver.advanceInTime();
				//---------------------


	            /*Teuchos::TimeMonitor solveTimeMonitor(*solveTime);

	            navierStokes->addBoundaries(bcFactory);
				navierStokes->addRhsFunction(rhs);						    
	            navierStokes->initializeProblem();
	            navierStokes->assemble();
	            navierStokes->setBoundariesRHS();

	            std::string nlSolverType = parameterListProblem->sublist("General").get("Linearization","FixedPoint");
	            NonLinearSolver<SC,LO,GO,NO> nlSolver( nlSolverType );
	            nlSolver.solve( *navierStokes );
	            comm->barrier();*/
	        }

			MAIN_TIMER_STOP(Solver);	


			MAIN_TIMER_START(Refinement," Step 3:	 meshRefinement");

			/* Refinement
			domainRefined.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );

			{

				ProblemPtr_Type problem = Teuchos::rcp_dynamic_cast<Problem_Type>( navierStokes , true);
				domainRefined = meshRefiner.globalAlgorithm( domainPressure,  domainVelocity, navierStokes->getSolution(), problem, rhs );
			}

			domainPressure = domainRefined;
			domainVelocity = domainPressure;
			
			j++;*/
			MAIN_TIMER_STOP(Refinement);	
        
        // ####################
       
            
      //  }

		MAIN_TIMER_STOP(Total);	
		Teuchos::TimeMonitor::report(cout,"Main");
   
    }

    return(EXIT_SUCCESS);
}



