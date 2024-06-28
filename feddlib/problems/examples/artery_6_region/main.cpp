#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/SCI.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_StackedTimer.hpp>

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;

void reactionTerm(double *, double *, double *);

int main(int argc, char *argv[])
{
    Teuchos::oblackhomestream blackhole;
    Teuchos::GlobalMPISession mpiSession;

    Teuchos::RCP<const Teuchos::Comm<int>> comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    Teuchos::CommandLineProcessor commandLineProcessor;
    string underlyingLibrary = "Tpetra";
    string simulationParametersXML = "simulationParameters.xml";
    string materialParametersXML = "materialParameters.xml";
    string solverParametersXML = "solverParameters.xml";
    string structurePreconditionerParametersXML = "preconditionerParameters_Structure.xml";
    string chemistryPreconditionerParametersXML = "preconditionerParameters_Chemistry.xml";
    commandLineProcessor.setOption("underlyingLibrary", &underlyingLibrary, "Underlying Library");
    commandLineProcessor.setOption("simulationsParameters", &simulationParametersXML, "xml file with simulation parameters");
    commandLineProcessor.setOption("materialParameters", &materialParametersXML, "xml file with material parameters");
    commandLineProcessor.setOption("solverParameters", &solverParametersXML, "xml file with solver parameters");
    commandLineProcessor.setOption("preconditionerParametersStructure", &structurePreconditionerParametersXML, "xml file with structure preconditoner parameters");
    commandLineProcessor.setOption("preconditionerParametersChemistry", &chemistryPreconditionerParametersXML, "xml file with chemistry preconditioner parameters");

    commandLineProcessor.recogniseAllOptions(true);
    commandLineProcessor.throwExceptions(true);

    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = commandLineProcessor.parse(argc, argv);
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    Teuchos::RCP<Teuchos::StackedTimer> stackedTimer = Teuchos::rcp(new StackedTimer("Structure-chemical interaction", true));
    bool verbose(comm->getRank() == 0);

    TimeMonitor::setStackedTimer(stackedTimer);
    {
        Teuchos::RCP<Teuchos::ParameterList> simulationParameters = Teuchos::getParametersFromXmlFile(simulationParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> materialParameters = Teuchos::getParametersFromXmlFile(materialParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> solverParameters = Teuchos::getParametersFromXmlFile(solverParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> structurePreconditionerParameters = Teuchos::getParametersFromXmlFile(structurePreconditionerParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> chemistryPreconditionersParamerters = Teuchos::getParanetersFromXmlFile(chemistryPreconditionerParametersXML);

        int dimension = simulationParameters->sublist("Simulation Parameters").get("Dimension");
        string discretizationType = simulationParameters->sublist("Simulation Parameters").get("Discretization");
        string preconditionerType = simulationParameters->sublist("Simulation Parameters").get("Preconditioner Type");

        Teuchos::RCP<Teuchos::ParameterList> allParameters = Teuchos::rcp(new Teuchos::ParameterList(*simulationParameters));

        Teuchos::RCP<Teuchos::ParameterList> preconditionerParameters = Teuchos::rcp(new Teuchos::ParameterList(*structurePreconditionerParameters));

        allParameters->setParameters(*materialParameters);

        Teuchos::RCP<Teuchos::ParameterList> allDiffusionParameters = Teuchos::rcp(new Teuchos::ParameterList(*chemistryPreconditionerParameters));
        allDiffusionParameters->sublist("Parameter")->setParameters(simulationParameters->sublist("Parameter Chem"));
        allDiffusionParameters->sublist("Parameter")->setParameters(simulationParameters->sublist("Simulation Parameters"));
        allDiffusionParameters->setParameters(*solverParameters);
        allDiffusionParameters->setParameters(*chemistryPreconditionerParameters);

        Teuchos::RCP<Teuchos::ParameterList> allStructureParameters = Teuchos::rcp(new Teuchos::ParameterList(*structurePreconditionerParameters));
        structurePreconditionerParameters->sublist("Parameter")->setParameters(simulationParameters->sublist("Parameter Solid"));

        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainP1Diffusion;
        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainP1Structure;
        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainP2Diffusion;
        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainP2Structure;
        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainDiffusion;
        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainStructure;

        domainP1Diffusion.reset(new FEDD::Domain<SC, LO, GO, NO>(comm, dimension));
        domainP1Structure.reset(new FEDD::Domain<SC, LO, GO, NO>(comm, dimension));
        domainP2Diffusion.reset(new FEDD::Domain<SC, LO, GO, NO>(comm, dimension));
        domainP2Structure.reset(new FEDD::Domain<SC, LO, GO, NO>(comm, dimension));

        std::vector<Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>>> domainP1Array(1);

        domainP1Array[0] = domainP1Structure;

        Teuchos::RCP<Teuchos::ParameterList> partitionerParameters = sublist(allParameters, "Mesh Partitioner");
        partitionerParameters->set("Build Edge List", true);
        partitionerParameters->set("Build Surface List", true);

        // Read the given mesh and partition it
        FEDD::MeshPartitioner<SC, LO, GO, NO> p1Partitioner(domainP1Array, partitionerParameters, "P1", dimension);
        p1Partitioner.readAndPartition(15);

        // Paraview export of node and element flags
        domainP1Structure->exportElementFlags();
        domainP1Structure->exportNodeFlags();

        // Build a P2 mesh using the P1 mesh
        domainDiffusion->buildP2ofP1Domain(domainP1Structure);
        domainStructure->buildP2ofP1Domain(domainP1Structure);

        // Set the degrees of freedem per node
        domainStructure->setDofs(dimension);
        domainDiffusion->setDofs(1);

        // Setting the reference configuration for both domains
        domainStructure->setReferenceConfiguration();
        domainDiffusion->setReferenceConfiguration();

        // Define the block system
        Teuchos::RCP<FEDD::SmallMatrix<int>> defTS;
        (*defTS)[0][0] = 1; // Structure
        (*defTS)[1][1] = 1; // Diffusion

        // Creating an object of the SCI system
        FEDD::SCI<SC, LO, GO, NO> sci(domainStructure, discretizationType,
                                      domainDiffusion, discretizationType,
                                      reactionTerm, allStructureParameters,
                                      allDiffusionParameters, allParameters,
                                      defTS);
        // Printing information
        sci.info();

        Teuchos::RCP<FEDD::BCBuilder<SC, LO, GO, NO>> bcFactory(new FEDD::BCBuilder<SC, LO, GO, NO>());
        Teuchos::RCP<FEDD::BCBuilder<SC, LO, GO, NO>> bcFactoryDiffusion(new FEDD::BCBuilder<SC, LO, GO, NO>());
        Teuchos::RCP<FEDD::BCBuilder<SC, LO, GO, NO>> bcFactoryStructure(new FEDD::BCBuilder<SC, LO, GO, NO>());

        // @TODO Set boundary conditions for plaque
    }
}

// @brief Reaction Term in the reaction diffusion equation
void reactionTerm(double *x, double *res, double *parameters)
{
    double m = 0.0;
    res[0] = m * x[0];
}