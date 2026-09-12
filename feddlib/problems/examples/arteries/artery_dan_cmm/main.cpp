// Comparison case for svMultiPhysics's Interface2/AceGen CCB element
// (tests/cases/def_diffu/artery_dan in the svMultiPhysics repo): the artery of
// ../paper_amlodipine/.../artery_dan (mesh, boundary conditions, pressure load)
// with the constrained-mixture (CMM) element (SCI_SMC_CMM_Active_Growth_Reorientation)
// in all seven tissue regions, set up so that both codes solve the same problem:
// the same P2 mesh (both build straight-edged P2 elements from the P1 file), the
// same Dirichlet node sets, the same material parameters and the same load/time
// stepping.
//
// svMultiPhysics prescribes Dirichlet conditions on faces, so the vertices that
// artery_dan pins (flags 13 and 14) become six-node pin faces there. The geometry
// override file (written by svMultiPhysics's tools/feddlib_mesh_to_svmp.py
// --dirichlet-override) sets the Dirichlet flags of exactly svMultiPhysics's
// constrained nodes:
//    2 - bottom face (z=0)                      -> Dirichlet_Z
//    3 - top face                               -> Dirichlet_Z
//    4 - outer wall                             -> concentration
//    5 - inner wall                             -> concentration (and the pressure, on its triangles)
//    6, 7 - bottom ring of the outer, inner wall -> Dirichlet_Z and concentration
//    8, 9 - top ring of the inner, outer wall    -> Dirichlet_Z and concentration
//   13 - pin face held in x and z               -> Dirichlet_X_Z
//   14 - pin faces held in y and z              -> Dirichlet_Y_Z
//   23, 24 - nodes of the pin faces 13, 14 on the outer wall -> as 13, 14, and concentration
//
// The pressure on the inner wall (flag 5) is ramped linearly as in
// svMultiPhysics's load.dat: lambda(t) = targetPressure/maxPressure * t / Ramp End Time,
// evaluated at t_{n+1} (artery_dan's ramp adds one load step: lambda(t + Load Step Size)).

#include "feddlib/core/General/BCBuilder.hpp"
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
#include <Teuchos_CommHelpers.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_StackedTimer.hpp>

#include <array>
#include <fstream>
#include <set>
#include <sstream>

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;

void reactionTerm(double *, double *, double *);
void loadFunction(double *, double *, double *);
void inflowChem(double *, double *, double, const double *);
void zeroDirichlet3D(double *, double *, double, const double *);

namespace {

struct GeometryOverrideEntry {
    std::array<double, 3> key;       // the node's straight-P2 position
    std::array<double, 3> position;  // new position (COORD entries)
    int flag;                        // new Dirichlet flag (FLAG entries)
    bool movesNode;
};

std::vector<GeometryOverrideEntry> readGeometryOverride(const std::string &fileName)
{
    std::ifstream in(fileName);
    TEUCHOS_TEST_FOR_EXCEPTION(!in, std::runtime_error, "Cannot open geometry override file " << fileName);
    std::vector<GeometryOverrideEntry> entries;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream ss(line);
        std::string kind;
        GeometryOverrideEntry e{};
        ss >> kind >> e.key[0] >> e.key[1] >> e.key[2];
        if (kind == "COORD")
        {
            e.movesNode = true;
            e.flag = -1;
            ss >> e.position[0] >> e.position[1] >> e.position[2];
        }
        else if (kind == "FLAG")
        {
            e.movesNode = false;
            ss >> e.flag;
        }
        else
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Unknown line in geometry override file: " << line);
        TEUCHOS_TEST_FOR_EXCEPTION(ss.fail(), std::runtime_error, "Malformed line in geometry override file: " << line);
        entries.push_back(e);
    }
    return entries;
}

// Matches entries to nodes by the node's current (straight P2) position, then moves
// nodes and sets flags. Nodes carrying one of the Dirichlet flags that are not listed
// lose it, so the constrained node sets are exactly the override's. Returns
// {#nodes moved, #flags set} for this point list.
std::array<int, 2> applyGeometryOverride(std::vector<std::vector<double>> &points, std::vector<int> &flags,
                                         const std::vector<GeometryOverrideEntry> &entries)
{
    const std::set<int> dirichletFlags = {2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 23, 24};
    const double tol2 = 1.e-16;
    std::vector<int> node(entries.size(), -1);
    for (size_t e = 0; e < entries.size(); e++)
        for (size_t i = 0; i < points.size(); i++)
        {
            double d2 = 0.;
            for (int c = 0; c < 3; c++)
            {
                const double d = points[i][c] - entries[e].key[c];
                d2 += d * d;
            }
            if (d2 < tol2)
            {
                node[e] = static_cast<int>(i);
                break;
            }
        }

    for (int &f : flags)
        if (dirichletFlags.count(f))
            f = 0;

    std::array<int, 2> applied = {0, 0};
    for (size_t e = 0; e < entries.size(); e++)
    {
        if (node[e] < 0)
            continue;
        if (entries[e].movesNode)
        {
            for (int c = 0; c < 3; c++)
                points[node[e]][c] = entries[e].position[c];
            applied[0]++;
        }
        else
        {
            flags[node[e]] = entries[e].flag;
            applied[1]++;
        }
    }
    return applied;
}

} // namespace

int main(int argc, char *argv[])
{
    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

    Teuchos::RCP<const Teuchos::Comm<int>> comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    Teuchos::CommandLineProcessor commandLineProcessor;
    std::string underlyingLibrary = "Tpetra";
    std::string simulationParametersXML = "simulationParameters.xml";
    std::string materialParametersXML = "materialParameters.xml";
    std::string solverParametersXML = "solverParameters.xml";
    std::string structurePreconditionerParametersXML = "preconditionerParameters_Structure.xml";
    std::string chemistryPreconditionerParametersXML = "preconditionerParameters_Chemistry.xml";
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

    Teuchos::RCP<Teuchos::StackedTimer> stackedTimer = Teuchos::rcp(new Teuchos::StackedTimer("Artery dan CMM comparison", true));
    bool verbose(comm->getRank() == 0);

    Teuchos::TimeMonitor::setStackedTimer(stackedTimer);
    {
        Teuchos::RCP<Teuchos::ParameterList> simulationParameters = Teuchos::getParametersFromXmlFile(simulationParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> materialParameters = Teuchos::getParametersFromXmlFile(materialParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> solverParameters = Teuchos::getParametersFromXmlFile(solverParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> structurePreconditionerParameters = Teuchos::getParametersFromXmlFile(structurePreconditionerParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> chemistryPreconditionerParamerters = Teuchos::getParametersFromXmlFile(chemistryPreconditionerParametersXML);

        int dimension = simulationParameters->sublist("Simulation Parameters").get("Dimension", 3);
        std::string discretizationType = simulationParameters->sublist("Simulation Parameters").get("Discretization", "P2");

        Teuchos::RCP<Teuchos::ParameterList> allParameters = Teuchos::rcp(new Teuchos::ParameterList(*simulationParameters));
        allParameters->sublist("Parameter").set("Chemistry Explicit", false);

        allParameters->setParameters(*materialParameters);
        allParameters->setParameters(*structurePreconditionerParameters);
        allParameters->setParameters(*solverParameters);

        Teuchos::RCP<Teuchos::ParameterList> allDiffusionParameters = Teuchos::rcp(new Teuchos::ParameterList(*chemistryPreconditionerParamerters));
        Teuchos::sublist(allDiffusionParameters, "Parameter")->setParameters(simulationParameters->sublist("Parameter Chem"));
        Teuchos::sublist(allDiffusionParameters, "Parameter")->setParameters(simulationParameters->sublist("Simulation Parameters"));
        allDiffusionParameters->setParameters(*solverParameters);
        allDiffusionParameters->setParameters(*chemistryPreconditionerParamerters);

        Teuchos::RCP<Teuchos::ParameterList> allStructureParameters = Teuchos::rcp(new Teuchos::ParameterList(*structurePreconditionerParameters));
        Teuchos::sublist(allStructureParameters, "Parameter")->setParameters(simulationParameters->sublist("Parameter Solid"));
        allStructureParameters->setParameters(*materialParameters);
        allStructureParameters->setParameters(*solverParameters);

        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainP1Structure;
        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainDiffusion;
        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainStructure;

        domainP1Structure.reset(new FEDD::Domain<SC, LO, GO, NO>(comm, dimension));
        domainDiffusion.reset(new FEDD::Domain<SC, LO, GO, NO>(comm, dimension));
        domainStructure.reset(new FEDD::Domain<SC, LO, GO, NO>(comm, dimension));

        std::vector<Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>>> domainP1Array(1);
        domainP1Array[0] = domainP1Structure;

        Teuchos::RCP<Teuchos::ParameterList> partitionerParameters = Teuchos::sublist(allParameters, "Mesh Partitioner");
        partitionerParameters->set("Build Edge List", true);
        partitionerParameters->set("Build Surface List", true);

        FEDD::MeshPartitioner<SC, LO, GO, NO> p1Partitioner(domainP1Array, partitionerParameters, "P1", dimension);
        p1Partitioner.readAndPartition(15);

        if (allParameters->sublist("General").get("ParaViewExport", true))
        {
            domainP1Structure->exportElementFlags();
            domainP1Structure->exportNodeFlags();
        }

        domainDiffusion->buildP2ofP1Domain(domainP1Structure);
        domainStructure->buildP2ofP1Domain(domainP1Structure);

        domainStructure->setDofs(dimension);
        domainDiffusion->setDofs(1);

        // Make the Dirichlet node sets identical to svMultiPhysics's before the
        // reference configuration is taken and the elements are built.
        std::string geometryOverride = partitionerParameters->get("Geometry Override", std::string(""));
        if (!geometryOverride.empty())
        {
            std::vector<GeometryOverrideEntry> entries = readGeometryOverride(geometryOverride);
            int expected[2] = {0, 0};
            for (const auto &e : entries)
                expected[e.movesNode ? 0 : 1]++;

            for (auto domain : {domainStructure, domainDiffusion})
            {
                applyGeometryOverride(*domain->getPointsRepeated(), *domain->getBCFlagRepeated(), entries);
                std::array<int, 2> local = applyGeometryOverride(*domain->getPointsUnique(), *domain->getBCFlagUnique(), entries);
                int global[2] = {0, 0};
                Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 2, local.data(), global);
                TEUCHOS_TEST_FOR_EXCEPTION(global[0] != expected[0] || global[1] != expected[1], std::runtime_error,
                                           "Geometry override " << geometryOverride << " matched " << global[0] << " of " << expected[0]
                                                                << " node moves and " << global[1] << " of " << expected[1] << " flags.");
            }
            if (verbose)
                std::cout << " Geometry override " << geometryOverride << ": set " << expected[1]
                          << " Dirichlet flags (svMultiPhysics's constrained nodes)." << std::endl;
        }

        domainStructure->setReferenceConfiguration();
        domainDiffusion->setReferenceConfiguration();

        Teuchos::RCP<FEDD::SmallMatrix<int>> defTS;
        defTS.reset(new FEDD::SmallMatrix<int>(2));
        (*defTS)[0][0] = 1;
        (*defTS)[1][1] = 1;

        std::vector<std::vector<double>> diffusionTensor(dimension, std::vector<double>(3));
        double D0 = allParameters->sublist("Parameter Diffusion").get("D0", 1.);
        for (int i = 0; i < dimension; i++)
        {
            diffusionTensor[0][0] = D0;
            diffusionTensor[1][1] = D0;
            diffusionTensor[2][2] = D0;
            if (i > 0)
            {
                diffusionTensor[i][i - 1] = 0;
                diffusionTensor[i - 1][i] = 0;
            }
            else
                diffusionTensor[i][i + 1] = 0;
        }

        FEDD::SCI<SC, LO, GO, NO> sci(domainStructure, discretizationType,
                                      domainDiffusion, discretizationType,
                                      diffusionTensor, reactionTerm,
                                      allStructureParameters, allDiffusionParameters,
                                      allParameters, defTS);
        sci.info();

        Teuchos::RCP<FEDD::BCBuilder<SC, LO, GO, NO>> bcFactory(new FEDD::BCBuilder<SC, LO, GO, NO>());
        Teuchos::RCP<FEDD::BCBuilder<SC, LO, GO, NO>> bcFactoryDiffusion(new FEDD::BCBuilder<SC, LO, GO, NO>());
        Teuchos::RCP<FEDD::BCBuilder<SC, LO, GO, NO>> bcFactoryStructure(new FEDD::BCBuilder<SC, LO, GO, NO>());

        double rampTimeStep = allParameters->sublist("Parameter").get("Load Step Size", 0.02);
        double timeRampEnd = allParameters->sublist("Parameter").get("Ramp End Time", 1.0);

        // Pressure in mmHg converted to kPa (1 mmHg = 0.133322 kPa), as in artery_dan.
        // The negative sign is FEDDLib's convention for an inflating internal pressure.
        double maxPressureMmHg = allParameters->sublist("Parameter").get("Max Pressure mmHg", 97.14);
        double targetPressureMmHg = allParameters->sublist("Parameter").get("Target Pressure mmHg", 85.0);
        double pressureReductionStartTime = allParameters->sublist("Parameter").get("Pressure Reduction Start Time", 1000000.0);
        double pressureReductionEndTime = allParameters->sublist("Parameter").get("Pressure Reduction End Time", 1000000.0);
        double pressureReductionAmountMmHg = allParameters->sublist("Parameter").get("Pressure Reduction Amount mmHg", 0.0);

        double pressure = -maxPressureMmHg * 0.133322;
        double initialLambda = targetPressureMmHg / maxPressureMmHg;
        double lambdaReduction = pressureReductionAmountMmHg / maxPressureMmHg;

        sci.problemStructureNonLin_->addParemeterRhs(pressure);
        sci.problemStructureNonLin_->addParemeterRhs(rampTimeStep);
        sci.problemStructureNonLin_->addParemeterRhs(timeRampEnd);
        sci.problemStructureNonLin_->addParemeterRhs(initialLambda);
        sci.problemStructureNonLin_->addParemeterRhs(pressureReductionStartTime);
        sci.problemStructureNonLin_->addParemeterRhs(pressureReductionEndTime);
        sci.problemStructureNonLin_->addParemeterRhs(lambdaReduction);
        sci.problemStructureNonLin_->addParemeterRhs( 0. ); // degree of the load function in space: develop's surface integral reads it from the last parameter

        sci.problemStructureNonLin_->addRhsFunction(loadFunction, 0);

        // Structure Dirichlet BCs (see the flag list at the top).
        std::vector<std::pair<int, std::string>> structureBCs = {
            {2, "Dirichlet_Z"}, {3, "Dirichlet_Z"}, {6, "Dirichlet_Z"}, {7, "Dirichlet_Z"},
            {8, "Dirichlet_Z"}, {9, "Dirichlet_Z"}, {13, "Dirichlet_X_Z"}, {23, "Dirichlet_X_Z"},
            {14, "Dirichlet_Y_Z"}, {24, "Dirichlet_Y_Z"}};
        for (const auto &bc : structureBCs)
        {
            bcFactoryStructure->addBC(zeroDirichlet3D, bc.first, 0, domainStructure, bc.second, dimension);
            bcFactory->addBC(zeroDirichlet3D, bc.first, 0, domainStructure, bc.second, dimension);
        }

        if (!sci.problemStructure_.is_null())
            sci.problemStructure_->addBoundaries(bcFactoryStructure);
        else
            sci.problemStructureNonLin_->addBoundaries(bcFactoryStructure);

        // Concentration on the inner (5, rings 7/8) and outer (4, rings 6/9, pin-face
        // nodes 23/24) walls: "Inflow Concentration" from "Inflow Start Time" on, zero before.
        std::vector<double> inflowParameters = {allParameters->sublist("Parameter").get("Inflow Start Time", 1.e7),
                                                allParameters->sublist("Parameter").get("Inflow Concentration", 2.0)};
        for (int wallFlag : {4, 5, 6, 7, 8, 9, 23, 24})
        {
            bcFactoryDiffusion->addBC(inflowChem, wallFlag, 0, domainDiffusion, "Dirichlet", 1, inflowParameters);
            bcFactory->addBC(inflowChem, wallFlag, 1, domainDiffusion, "Dirichlet", 1, inflowParameters);
        }
        sci.problemChem_->addBoundaries(bcFactoryDiffusion);

        sci.addBoundaries(bcFactory);

        sci.initializeProblem();
        sci.initializeCE();

        FEDD::DAESolverInTime<SC, LO, GO, NO> daeTimeSolver(allParameters, comm);
        daeTimeSolver.defineTimeStepping(*defTS);
        daeTimeSolver.setProblem(sci);
        sci.assemble();
        daeTimeSolver.setupTimeStepping();
        daeTimeSolver.advanceInTime();
    }
    FEDD::TimeMonitor_Type::report(std::cout);
    stackedTimer->stop("Artery dan CMM comparison");
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stackedTimer->report((std::cout), comm, options);

    return (EXIT_SUCCESS);
}

// No separate reaction term outside the AceGen kernel's own internal Rc computation
// (m=0 is a no-op) -- matches svMultiPhysics's assembly.
void reactionTerm(double *x, double *res, double *parameters)
{
    double m = 0.0;
    res[0] = m * x[0];
}

/* Parameter order (see Problem_decl.hpp addParemeterRhs / FE_def.hpp
 * assemblySurfaceIntegralExternal):
 * parameters[0]: current time (t_{n+1})
 * parameters[1]: pressure (kPa)
 * parameters[2]: rampTimeStep (not needed by the ramp below)
 * parameters[3]: timeRampEnd
 * parameters[4]: initialLambda
 * parameters[5]: pressureReductionStartTime
 * parameters[6]: pressureReductionEndTime
 * parameters[7]: lambdaReduction
 * parameters[8]: degree of the load function in space (0)
 * parameters[9]: surface flag
 */
void loadFunction(double *x, double *res, double *parameters)
{
    res[0] = 0.0;

    double currentTime = parameters[0];
    double pressure = parameters[1];
    double timeRampEnd = parameters[3];
    double initialLambda = parameters[4];
    double pressureReductionStartTime = parameters[5];
    double pressureReductionEndTime = parameters[6];
    double lambdaReduction = parameters[7];
    double surfaceFlag = parameters[9];

    double lambda = 0.0;
    double currentLambdaReduction = 0.0;

    // currentTime is already t_{n+1} (DAESolverInTime advances time before the solve),
    // so lambda(t) = initialLambda * t / timeRampEnd gives svMultiPhysics's load.dat ramp.
    if (currentTime < timeRampEnd)
        lambda = initialLambda * currentTime / timeRampEnd;
    else
        lambda = initialLambda;

    if (currentTime >= pressureReductionStartTime && currentTime <= pressureReductionEndTime) {
        double reductionProgress = (currentTime - pressureReductionStartTime) /
                                 (pressureReductionEndTime - pressureReductionStartTime);
        currentLambdaReduction = lambdaReduction * reductionProgress;
    } else if (currentTime > pressureReductionEndTime) {
        currentLambdaReduction = lambdaReduction;
    }

    lambda = lambda - currentLambdaReduction;
    if (lambda < 0.0) lambda = 0.0;

    if (surfaceFlag == 5) // inner wall
        res[0] = pressure * lambda;
}

// parameters[0]: start time, parameters[1]: concentration from then on.
void inflowChem(double *x, double *res, double t, const double *parameters)
{
    res[0] = (t >= parameters[0]) ? parameters[1] : 0.;
}

// Fix all degrees of freedom (BCBuilder picks out only the component(s) its BC type
// string names, e.g. "Dirichlet_Z" only uses res[2]).
void zeroDirichlet3D(double *x, double *res, double t, const double *parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
}
