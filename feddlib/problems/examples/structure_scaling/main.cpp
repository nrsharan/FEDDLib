#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/specific/NonLinElasAssFE.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_StackedTimer.hpp>

void rhsDummy2D(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    return;
}

void rhsDummy(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void zeroBC(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
}

void zeroDirichlet2D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;

    return;
}


void zeroDirichlet(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
}

void reactionFunc(double* x, double* res, double* parameters){
	
    double m = 0.0;	
    res[0] = m * x[0];

}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void inflowChem(double* x, double* res, double t, const double* parameters)
{
	if(t>=parameters[0])
    	res[0] = 1.;
    else	
    	res[0] = 0.;
    return;
}


void rhsX(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    
    res[0] = parameters[1];
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void rhsY(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] =  parameters[1];
    res[2] = 0.;
    return;
}

void rhsZ(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    res[2] = parameters[1];
    return;
}

// Parameter Structure
// 0 : time
// 1 : force
// 2 : loadStep (lambda)
// 3 : LoadStep end time
// 4 : Flag 

void rhsYZ(double* x, double* res, double* parameters){

    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];

  	res[0] =0.;
    res[1] =0.;
    res[2] =0.;

    if(parameters[0]+1.e-12 < TRamp)
        force = (parameters[0]+loadStepSize) * parameters[1] / TRamp ;
    else
        force = parameters[1];


    if(parameters[5] == 4  || parameters[5] == 5){
      	res[0] = force;
        res[1] = force;
        res[2] = force;
    }
    
    return;
}


// Parameter Structure
// 0 : time
// 1 : force
// 2 : loadStep (lambda)
// 3 : LoadStep end time
// 4 : Flag 

void rhsArtery(double* x, double* res, double* parameters){

    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];

  	res[0] =0.;
    res[1] =0.;
    res[2] =0.;
    
    if(parameters[0]+1.e-12 < TRamp)
        force = (parameters[0]+loadStepSize) * parameters[1] / TRamp ;
    else
        force = parameters[1];


    if(parameters[5] == 5){
      	res[0] = force;
        res[1] = force;
        res[2] = force;
    }
    
    return;
}

void dummyFunc(double* x, double* res, double t, const double* parameters)
{
    return;
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
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
   
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");    
    
    string xmlProblemStructureFile = "parametersProblemStructure.xml";  
    myCLP.setOption("problemfileStructure",&xmlProblemStructureFile,".xml file with Inputparameters.");    
 
    string xmlSolverFile = "parametersSolver.xml"; 
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");
 
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
        mpiSession.~GlobalMPISession();
        return 0;
    }
	Teuchos::RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Structure-chemical interaction",true));
    bool verbose (comm->getRank() == 0);
    TimeMonitor::setStackedTimer(stackedTimer);
    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
       
        ParameterListPtr_Type parameterListProblemStructure = Teuchos::getParametersFromXmlFile(xmlProblemStructureFile);
        
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",2);
        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","unstructured");
        int 		m				= parameterListProblem->sublist("Parameter").get("H/h",5);
        string      FEType        = parameterListProblem->sublist("Parameter").get("Discretization","P2");
        string precMethod = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        int         n;

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;   
          
		ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
       

        parameterListAll->setParameters(*parameterListSolver);
        parameterListAll->setParameters(*parameterListPrec);
                    
    	parameterListAll->setParameters(*parameterListProblemStructure);
              
        ParameterListPtr_Type parameterListStructureAll(new Teuchos::ParameterList(*parameterListPrec));
        sublist(parameterListStructureAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Solid") );
        parameterListStructureAll->setParameters(*parameterListPrec);
        parameterListStructureAll->setParameters(*parameterListProblem);
        parameterListStructureAll->setParameters(*parameterListProblemStructure);
		
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
            cout << "############ Starting SCI  ... ################" <<endl;
            cout << "###############################################" <<endl;
        }


        DomainPtr_Type domainP1struct;
        DomainPtr_Type domainP2struct;
        
        DomainPtr_Type domainStructure;
        
        std::string bcType = parameterListAll->sublist("Parameter").get("BC Type","Cube");
        
        std::string rhsType = parameterListAll->sublist("Parameter").get("RHS Type","Constant");
    
        int minNumberSubdomains=1;
       
        if (!meshType.compare("structured")) {
		    TEUCHOS_TEST_FOR_EXCEPTION( size%minNumberSubdomains != 0 , std::logic_error, "Wrong number of processors for structured mesh.");
	  	    n = (int)(std::pow( size/minNumberSubdomains, 1/3.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
	        std::vector<double> x(3);
	        x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
	        domainStructure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));

		    domainStructure->buildMesh( 3,"Square", dim, FEType, n, m, numProcsCoarseSolve);
		}
        else if (!meshType.compare("unstructured")) {
        

		    domainP1struct.reset( new Domain_Type( comm, dim ) );
		    domainP2struct.reset( new Domain_Type( comm, dim ) );
            
            MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
            domainP1Array[0] = domainP1struct;
            
            ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );                    
            
            pListPartitioner->set("Build Edge List",true);
		    pListPartitioner->set("Build Surface List",true);
		                    
		    MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
		    
		    int volumeID=10;
		    if(bcType=="Artery" || bcType == "Artery Full" || bcType == "Artery Plaque")
		    	volumeID = 15;
		    else if(bcType=="Realistic Artery 1" || bcType=="Realistic Artery 2")
		    	volumeID = 21;
		    	
		    partitionerP1.readAndPartition(volumeID);
		    
		    
            if (!FEType.compare("P2")){
				domainP2struct->buildP2ofP1Domain( domainP1struct );
				domainStructure = domainP2struct;   
			}        
			else{
                TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only P2 discretization allowed");                               

				domainStructure = domainP1struct;
			}
        }
        domainStructure->setDofs(dim);
     
        if (parameterListAll->sublist("General").get("ParaView export subdomains",false) ){
		   // ########################
		    // Flags check
		    // ########################

			Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaF(new ExporterParaView<SC,LO,GO,NO>());

			Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution(new MultiVector<SC,LO,GO,NO>(domainStructure->getMapUnique()));
			vec_int_ptr_Type BCFlags = domainStructure->getBCFlagUnique();

			Teuchos::ArrayRCP< SC > entries  = exportSolution->getDataNonConst(0);
			for(int i=0; i< entries.size(); i++){
				entries[i] = BCFlags->at(i);
			}

			Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst = exportSolution;

			exParaF->setup("Flags", domainStructure->getMesh(), FEType);

			exParaF->addVariable(exportSolutionConst, "Flags", "Scalar", 1,domainStructure->getMapUnique(), domainStructure->getMapUniqueP2());

			exParaF->save(0.0);


            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaE(new ExporterParaView<SC,LO,GO,NO>());

			Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolutionE(new MultiVector<SC,LO,GO,NO>(domainStructure->getElementMap()));
			
			Teuchos::ArrayRCP< SC > entriesE  = exportSolutionE->getDataNonConst(0);
			for(int i=0; i<domainStructure->getElementsC()->numberElements(); i++){
				entriesE[i] = domainStructure->getElementsC()->getElement(i).getFlag();
			}

			Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConstE = exportSolutionE;

			exParaE->setup("Flags_Elements", domainStructure->getMesh(), "P0");

			exParaE->addVariable(exportSolutionConstE, "Flags_Elements", "Scalar", 1,domainStructure->getElementMap());

			exParaE->save(0.0);
		
	

            
            if (verbose)
                std::cout << "\t### Exporting subdomains ###\n";

            typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
            typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
            typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
            typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
            typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
            // Same subdomain for solid and chemistry, as they have same domain
            {
               MultiVectorPtr_Type vecDecomposition = rcp(new MultiVector_Type( domainStructure->getElementMap() ) );
                MultiVectorConstPtr_Type vecDecompositionConst = vecDecomposition;
                vecDecomposition->putScalar(comm->getRank()+1.);
                
                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
                
                exPara->setup( "subdomains_solid", domainStructure->getMesh(), "P0" );
                
                exPara->addVariable( vecDecompositionConst, "subdomains", "Scalar", 1, domainStructure->getElementMap());
                exPara->save(0.0);
                exPara->closeExporter();
            }
           
        
     
        }

       
        domainStructure->setReferenceConfiguration();              
			

        NonLinElasAssFE<SC,LO,GO,NO> nonLinElas(domainStructure,FEType,parameterListAll);
        
        nonLinElas.info();
        
            
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) ); 
                   
    
        if(dim == 2)
        {
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only 3D Test available");                               
        }
        else if(dim == 3 && bcType=="Cube")
        {

            bcFactory->addBC(zeroDirichlet3D, 1, 0, domainStructure, "Dirichlet_X", dim);
            bcFactory->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Y", dim);
            bcFactory->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
            
            bcFactory->addBC(zeroDirichlet3D, 0, 0, domainStructure, "Dirichlet", dim);
            bcFactory->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_X_Y", dim);
            bcFactory->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Y_Z", dim);
            bcFactory->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_X_Z", dim);
            
        }
        else if(dim==3 && bcType=="Artery"){
        
			bcFactory->addBC(zeroDirichlet3D, 1, 0, domainStructure, "Dirichlet_Y", dim);
			bcFactory->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_X", dim);
			bcFactory->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 4, 0, domainStructure, "Dirichlet_Z", dim);
			
			bcFactory->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Z", dim);

			bcFactory->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Y_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_X_Z", dim);

			bcFactory->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_X", dim);
			bcFactory->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet_Y", dim);

			bcFactory->addBC(zeroDirichlet3D, 11, 0, domainStructure, "Dirichlet_Y_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 12, 0, domainStructure, "Dirichlet_X_Z", dim);
    
        }
       
        

       
        // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
            
             
        if(bcType=="Cube"){
			if(rhsType=="Constant")
    		 	nonLinElas.addRhsFunction( rhsYZ,0 );
    		/*f(rhsType=="Paper")
        		nonLinElas.problemStructureNonLin_->addRhsFunction( rhsCubePaper,0 );
    		if(rhsType=="Heart Beat")
        		nonLinElas.problemStructureNonLin_->addRhsFunction( rhsHeartBeatCube,0 );*/
			
		}
		else if(bcType=="Artery"){
			if(rhsType=="Constant")
    		 	nonLinElas.addRhsFunction( rhsArtery,0 );
			/*if(rhsType=="Paper")
    		 	nonLinElas.problemStructureNonLin_->addRhsFunction( rhsArteryPaper,0 );
    		if(rhsType=="Heart Beat")
    		 	nonLinElas.problemStructureNonLin_->addRhsFunction( rhsHeartBeatArtery,0 );
            if(rhsType=="Paper Pulse")
    		 	nonLinElas.problemStructureNonLin_->addRhsFunction(rhsArteryPaperPulse,0 );
            if(rhsType=="Heart Beat Pulse")
    		 	nonLinElas.problemStructureNonLin_->addRhsFunction(rhsHeartBeatArteryPulse,0 );*/
		}
        
        double force = parameterListAll->sublist("Parameter").get("Volume force",1.);
    
        nonLinElas.addParemeterRhs( force );
        double loadStep = parameterListAll->sublist("Parameter").get("Load Step Size",1.);
        double loadRampEnd= parameterListAll->sublist("Parameter").get("Load Ramp End",1.);
        nonLinElas.addParemeterRhs( loadStep );
        nonLinElas.addParemeterRhs( loadRampEnd );
        double heartBeatStart= parameterListAll->sublist("Parameter").get("Heart Beat Start",70.);
        nonLinElas.addParemeterRhs( heartBeatStart );

                  
        // #####################
        // Problem
        // #####################
        nonLinElas.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

        nonLinElas.initializeProblem();

        // Matrizen assemblieren
        nonLinElas.assemble();
                    
                    
      	// ######################
        // Zeitintegration
        // ######################
        DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);

        // Only one block for structural problem
        SmallMatrix<int> defTS(1);
        defTS[0][0] = 1;

        // Uebergebe auf welchen Bloecken die Zeitintegration durchgefuehrt werden soll
        // und Uebergabe der parameterList, wo die Parameter fuer die Zeitintegration drin stehen
        daeTimeSolver.defineTimeStepping(defTS);

        // Uebergebe das (nicht) lineare Problem
        daeTimeSolver.setProblem(nonLinElas);

        // Setup fuer die Zeitintegration, wie z.B. Aufstellen der Massematrizen auf den Zeilen, welche in
        // defTS definiert worden sind.
        daeTimeSolver.setupTimeStepping();

        // Fuehre die komplette Zeitintegration + Newton + Loesen + Exporter durch
        daeTimeSolver.advanceInTime();
    }
    TimeMonitor_Type::report(std::cout);
        stackedTimer->stop("Structure-chemical interaction");
	StackedTimer::OutputOptions options;
	options.output_fraction = options.output_histogram = options.output_minmax = true;
	stackedTimer->report((std::cout),comm,options);
	
    return(EXIT_SUCCESS);
}
  
