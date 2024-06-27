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

        Teuchos::RCP<Teuchos::ParameterList> allChemistryParameters = Teuchos::rcp(new Teuchos::ParameterList(*chemistryPreconditionerParameters));
        allChemistryParameters->sublist("Parameter")->setParameters(simulationParameters->sublist("Parameter Chem"));
        allChemistryParameters->sublist("Parameter")->setParameters(simulationParameters->sublist("Simulation Parameters"));
        allChemistryParameters->setParameters(*solverParameters);
        allChemistryParameters->setParameters(*chemistryPreconditionerParameters);

        Teuchos::RCP<Teuchos::ParameterList> allStructureParameters = Teuchos::rcp(new Teuchos::ParameterList(*structurePreconditionerParameters));
        structurePreconditionerParameters->sublist("Parameter")->setParameters(simulationParameters->sublist("Parameter Solid"));

        Teuchos::RCP<FEDDLib::Domain<SC,LO,GO,NO>> domainP1Diffusion;
        Teuchos::RCP<FEDDLib::Domain<SC,LO,GO,NO>> domainP1Structure;
        Teuchos::RCP<FEDDLib::Domain<SC,LO,GO,NO>> domainP2Diffusion;
        Teuchos::RCP<FEDDLib::Domain<SC,LO,GO,NO>> domainP2Structure;

        Teuchos::RCP<FEDDLib::Domain<SC,LO,GO,NO>> domainDiffusion;
        Teuchos::RCP<FEDDLib::Domain<SC,LO,GO,NO>> domainStructure;




        

    }
}
