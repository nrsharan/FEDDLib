#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/LinElas.hpp"
#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include "feddlib/problems/specific/LinElasAssFE.hpp"
#include "feddlib/problems/specific/NonLinElasAssFE.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include "feddlib/problems/Solver/NonLinearSolver.hpp"

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

	//cout << " Force: " << force << " flag = " << parameters[5] << endl;
    if(parameters[5] == 5){
      	res[0] = force;
        res[1] = force;
        res[2] = force;
    }
    
    return;
}

// Parameter Structure
// 0 : time
// 1 : force
// 2 : loadStepSize
// 3 : LoadStep end time
// 4 : Flag 

void rhsHeartBeatCube(double* x, double* res, double* parameters){

    res[0] =0.;
    res[1] =0.;
    res[2] = 0.;
    double lambda=0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double heartBeatStart = parameters[4];
    
	double a0    = 11.693284502463376;
	double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
		  0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
		  0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
	double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
		  -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
		  0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
		         
    double Q = 0.5*a0;
    

    double t_min = parameters[0] - fmod(parameters[0],1.0); //FlowConditions::t_start_unsteady;
    double t_max = t_min + 0.52; // One heartbeat lasts 1.0 second    
    double y = M_PI * ( 2.0*( parameters[0]-t_min ) / ( t_max - t_min ) -0.87  );
    
    for(int i=0; i< 20; i++)
        Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
    
    
    // Remove initial offset due to FFT
    Q -= 0.026039341343493;
    Q = (Q - 2.85489)/(7.96908-2.85489);
    
    bool Qtrue=false;
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0]+1.0e-10 < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))+1.e-10< 0.5)
    	lambda= 0.75;
    else{
        lambda = 0.75+0.25*Q;//*0.005329; // 0.775+0.125 * cos(4*M_PI*(parameters[0]));
        Qtrue = true; 
    } 
  
    double forceDirection = force/fabs(force);
    if(parameters[5]==5 || parameters[5]==4){
        res[0] =lambda*force;//+forceDirection*Q;
        res[1] =lambda*force;//+forceDirection*Q;
        res[2] =lambda*force;//+forceDirection*Q;        
    } 
    
   /* if(parameters[0]< heartBeatStart){
    	Q = 0.;
    }
    
    if(parameters[0]+1e-12 < TRamp)
        force = force * (parameters[0]+loadStepSize);
    
    if(parameters[5] == 5 || parameters[5] == 4){
     	res[0] = force+Q*0.005329;
        res[1] = force+Q*0.005329;
       	res[2] = force+Q*0.005329;
       	      	
    }*/
      
}
// Parameter Structure
// 0 : time
// 1 : force
// 2 : loadStepSize
// 3 : LoadStep end time
// 4 : Flag 
void rhsHeartBeatArtery(double* x, double* res, double* parameters){

    res[0] =0.;
    res[1] =0.;
    res[2] = 0.;
    double lambda=0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double heartBeatStart = parameters[4];
    
	double a0    = 11.693284502463376;
	double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
		  0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
		  0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
	double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
		  -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
		  0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
		         
    double Q = 0.5*a0;
    

    double t_min = parameters[0] - fmod(parameters[0],1.0); //FlowConditions::t_start_unsteady;
    double t_max = t_min + 0.5; // One heartbeat lasts 1.0 second    
    double y = M_PI * ( 2.0*( parameters[0]-t_min ) / ( t_max - t_min )-1.);// -0.87  );
    
    for(int i=0; i< 20; i++)
        Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
    
    
    // Remove initial offset due to FFT
    Q -= 0.026039341343493;
    Q = (Q - 2.85489)/(7.96908-2.85489);
    
    bool Qtrue=false;
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0]+1.0e-10 < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))+1.e-10< 0.5)
    	lambda= 0.75;
    else{
        lambda = 0.75+0.25*Q;//*0.005329; // 0.775+0.125 * cos(4*M_PI*(parameters[0]));
        Qtrue = true; 
    } 
  
    double forceDirection = force/fabs(force);
    if(parameters[5]==5){
        res[0] =lambda*force;//+forceDirection*Q;
        res[1] =lambda*force;//+forceDirection*Q;
        res[2] =lambda*force;//+forceDirection*Q;        
    } 
    
   /* if(parameters[0]< heartBeatStart){
    	Q = 0.;
    }
    
    if(parameters[0]+1e-12 < TRamp)
        force = force * (parameters[0]+loadStepSize);
    
    if(parameters[5] == 5 || parameters[5] == 4){
     	res[0] = force+Q*0.005329;
        res[1] = force+Q*0.005329;
       	res[2] = force+Q*0.005329;
       	      	
    }*/
      
}

void rhsHeartBeatArteryPulse(double* x, double* res, double* parameters){

    res[0] =0.;
    res[1] =0.;
    res[2] = 0.;
    double lambda=0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double heartBeatStart = parameters[4];
    
	    
    bool Qtrue=false;
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0]+1.0e-10 < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0]+1.0e-10 < heartBeatStart + 1.0)
		lambda = 0.75;
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))+1.e-10> 0.6)
    	lambda= 0.75;
    else{ // ( parameters[0] >= heartBeatStart + 0.5){ [0,0.6]
        // Within one second the heart beat passes through the artery
        double t= parameters[0] - std::floor(parameters[0]); // [0.5,1] Intervall
        double z= x[2];
        double t_z = t-1./80.*z; //  x  = 0.1 / 8 = 0.0125

        if(t_z > 0. && t_z < 0.5){
            double a0    = 11.693284502463376;
            double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
                0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
                0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
            double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
                -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
                0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
                        
            double Q = 0.5*a0;
        

            double t_min = 0; //parameters[0] - fmod(parameters[0],1.0); //FlowConditions::t_start_unsteady;
            double t_max = t_min + 0.5; // One heartbeat lasts 1.0 second    
            double y = M_PI * ( 2 * ( t_z-t_min ) / ( t_max - t_min ) - 1 );
           
            
            for(int i=0; i< 20; i++)
                Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
            
            
            // Remove initial offset due to FFT
            Q -= 0.026039341343493;
            Q = (Q - 2.85489)/(7.96908-2.85489);
            if(Q < 0 )
                Q = 0.;

            lambda = 0.75+0.25*Q;

        }
        else
            lambda = 0.75;
        
        if(lambda<0.75+1.e-12)
            lambda=0.75;

    } 
  
    double forceDirection = force/fabs(force);
    if(parameters[5]==5){
        res[0] =lambda*force;//+forceDirection*Q;
        res[1] =lambda*force;//+forceDirection*Q;
        res[2] =lambda*force;//+forceDirection*Q;        
    } 
    
}

// Parameter Structure
// 0 : time
// 1 : force
// 2 : loadStepSize
// 3 : LoadStep end time
// 4 : Flag 
void rhsArteryPaperPulse(double* x, double* res, double* parameters){

    res[0] =0.;
    res[1] =0.;
    res[2] = 0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double lambda=0.;
    double heartBeatStart = parameters[4];
    
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0] < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))< 0.5)
    	lambda= 0.75;
    else{
        double tinc = parameters[0] - std::floor(parameters[0]);
        double Q = -sin(1/16.*M_PI*x[2]-M_PI*(tinc-0.5)*3.0);
        if(Q < 0.+1.e-12){
            Q = 0.;
            lambda=0.75;
        }
        else{
            lambda =0.75+0.25*Q;//0.875 - 0.125
        }
    }
    

    if(parameters[5]==5){
        res[0] =lambda*force;
        res[1] =lambda*force;
        res[2] =lambda*force; 
        
       if(fabs(lambda*force)<0.75*0.016)
        cout << " ALARMAAAAA lamba=" << lambda << " force=" << force << " t= " << parameters[0] << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ----------------" << endl;
    }
      
}

void rhsArteryPaper(double* x, double* res, double* parameters){

    res[0] =0.;
    res[1] =0.;
    res[2] = 0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double lambda=0.;
    double heartBeatStart = parameters[4];
    
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0] < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))< 0.5)
    	lambda= 0.75;
    else
        lambda = 0.875 - 0.125 * cos(4*M_PI*(parameters[0]));
     
 
    if(parameters[5]==5){
        res[0] =lambda*force;
        res[1] =lambda*force;
        res[2] =lambda*force; 
        
       
    }
      
}

void rhsCubePaper(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[2] = 0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double lambda=0.;
    
    double heartBeatStart = parameters[4];
    
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0] < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))< 0.5)
    	lambda= 0.75;
    else
        lambda = 0.875 - 0.125 * cos(4*M_PI*(parameters[0]));
     
     
    if(parameters[5] == 5 || parameters[5] == 4){
        res[0] =lambda*force;
        res[1] =lambda*force;
        res[2] =lambda*force; 
    }
    
    

            
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
     
    string xmlSolverFile = "parametersSolver.xml"; 
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");
    
    string xmlPrecFileStructure = "parametersPrec.xml";
    myCLP.setOption("precfileStructure",&xmlPrecFileStructure,".xml file with Inputparameters.");
    string xmlPrecFileChem = "parametersPrecChem.xml";
    myCLP.setOption("precfileChem",&xmlPrecFileChem,".xml file with Inputparameters.");
    
 	string xmlBlockPrecFile = "parametersPrecBlock.xml";
    myCLP.setOption("blockprecfile",&xmlBlockPrecFile,".xml file with Inputparameters.");
   
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
    bool verbose (comm->getRank() == 0);

    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
       
        ParameterListPtr_Type parameterListProblemStructure = Teuchos::getParametersFromXmlFile(xmlProblemStructureFile);
        
        ParameterListPtr_Type parameterListSolver= Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListPrecStructure = Teuchos::getParametersFromXmlFile(xmlPrecFileStructure);
  
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);

 		int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",2);
        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","unstructured");
        int 		m				= parameterListProblem->sublist("Parameter").get("H/h",5);
        string      discType        = parameterListProblem->sublist("Parameter").get("Discretization","P2");
        string precMethod = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        int         n;
       
        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;     
        
        parameterListAll->setParameters(*parameterListSolver);

        parameterListAll->setParameters(*parameterListPrec);
        parameterListAll->setParameters(*parameterListProblemStructure);
        
        
        ParameterListPtr_Type parameterListStructureAll(new Teuchos::ParameterList(*parameterListPrecStructure));
        sublist(parameterListStructureAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Solid") );
        parameterListStructureAll->setParameters(*parameterListPrecStructure);
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
            cout << "############ Starting  ... ################" <<endl;
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
		    /*if (dim == 2) {
		        n = (int) (std::pow( size/minNumberSubdomains ,1/2.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
		        std::vector<double> x(2);
		        x[0]=0.0;    x[1]=0.0;
		        domainStructure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
		        domainChem.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
		    }
		    else if (dim == 3){*/
		        n = (int)(std::pow( size/minNumberSubdomains, 1/3.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
		        std::vector<double> x(3);
		        x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
		        domainStructure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
		    //}
		    domainStructure->buildMesh( 3,"Square", dim, discType, n, m, numProcsCoarseSolve);
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
		    if(bcType=="Artery")
		    	volumeID = 15;
		    
		    partitionerP1.readAndPartition(volumeID);
		    
		    
            if (!discType.compare("P2")){
				//domainP2chem->buildP2ofP1Domain( domainP1struct );
				domainP2struct->buildP2ofP1Domain( domainP1struct );
				
				domainStructure = domainP2struct;   
			}        
			else{
				domainStructure = domainP1struct;
			}
        }
        domainStructure->setDofs(dim);
  
       
        domainStructure->setReferenceConfiguration();
       
        Teuchos::RCP<SmallMatrix<int>> defTS;

        defTS.reset( new SmallMatrix<int> (1) );

        // Stucture
        (*defTS)[0][0] = 0.;
       
        NonLinElasAssFE<SC,LO,GO,NO> NonLinElasAssFE( domainStructure, discType, parameterListAll );      
            
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) ); 
                  
    
        // Struktur-RW
        

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
        else if(dim==3 && bcType=="Artery Full"){
        
			bcFactory->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
            bcFactory->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Z", dim);
			
			bcFactory->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_X_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Y_Z", dim);

			

		
        }
		
       
        // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
        NonLinElasAssFE.addBoundaries(bcFactory);
        
        // RHS dummy for structure
            
                
        if(bcType=="Cube"){
			if(rhsType=="Constant")
    		 	NonLinElasAssFE.addRhsFunction( rhsYZ,0 );
    		if(rhsType=="Paper")
        		NonLinElasAssFE.addRhsFunction( rhsCubePaper,0 );
    		if(rhsType=="Heart Beat")
        		NonLinElasAssFE.addRhsFunction( rhsHeartBeatCube,0 );
			
		}
		else if(bcType=="Artery" || bcType == "Artery Full" ){
			if(rhsType=="Constant")
    		 	NonLinElasAssFE.addRhsFunction( rhsArtery,0);
			if(rhsType=="Paper")
    		 	NonLinElasAssFE.addRhsFunction( rhsArteryPaper,0 );
    		if(rhsType=="Heart Beat")
    		 	NonLinElasAssFE.addRhsFunction( rhsHeartBeatArtery,0 );
            if(rhsType=="Paper Pulse")
    		 	NonLinElasAssFE.addRhsFunction(rhsArteryPaperPulse,0 );
            if(rhsType=="Heart Beat Pulse")
    		 	NonLinElasAssFE.addRhsFunction(rhsHeartBeatArteryPulse,0 );
		}
        
        double force = parameterListAll->sublist("Parameter").get("Volume force",1.);
        if(bcType == "Artery Full")
            force = force * -1.;
        NonLinElasAssFE.addParemeterRhs( force );
        double loadStep = parameterListAll->sublist("Parameter").get("Load Step Size",1.);
        double loadRampEnd= parameterListAll->sublist("Parameter").get("Load Ramp End",1.);
        NonLinElasAssFE.addParemeterRhs( loadStep );
        NonLinElasAssFE.addParemeterRhs( loadRampEnd );
        double heartBeatStart= parameterListAll->sublist("Parameter").get("Heart Beat Start",70.);
        NonLinElasAssFE.addParemeterRhs( heartBeatStart );

           
                 

        // #####################
        // Zeitintegration
        // #####################
        NonLinElasAssFE.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

        NonLinElasAssFE.initializeProblem();
        // Matrizen assemblieren
        NonLinElasAssFE.assemble();
                    

                    
        DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);

        // Uebergebe auf welchen Bloecken die Zeitintegration durchgefuehrt werden soll
        // und Uebergabe der parameterList, wo die Parameter fuer die Zeitintegration drin stehen
        daeTimeSolver.defineTimeStepping(*defTS);

        // Uebergebe das (nicht) lineare Problem
        daeTimeSolver.setProblem(NonLinElasAssFE);

        // Setup fuer die Zeitintegration, wie z.B. Aufstellen der Massmatrizen auf den Zeilen, welche in
        // defTS definiert worden sind.
        daeTimeSolver.setupTimeStepping();

        daeTimeSolver.advanceInTime();
    }
    	    TimeMonitor_Type::report(std::cout);

    return(EXIT_SUCCESS);
}

