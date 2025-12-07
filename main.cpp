#include <climits> // For INT_MAX
#include <ctime>   // For CLOCKS_PER_SEC
#include <set>
#include "main.h"

int main(int argc, const char *argv[]) {
    // Set parameters from stdin.
    Params params = getParams(argc, argv);

    // Read the problem instance from input file.
    SVRPInstance instance;
    if (!filehandler::readInstance(params, instance)) {
        cerr << "Could not read input file " << params.inputFile << endl;
        std::exit(1);
    }
    instance.print();

    // Run the algorithm.
    auto started = chrono::high_resolution_clock::now();

    // If using In-Out separation, we first run the base algorithm.
    Params newParams = params;
    SVRPSolution heuristicSolution;
    if (params.inOut) {
        newParams.timeLimit = 30;
        newParams.paradaSet = false;
        newParams.sriCuts = true;
        newParams.basicVRPSD = false;
        newParams.scenarioOptimalPRCuts = false;
        NodeRecourseModel nodeRecourseModel1(instance, newParams);
        nodeRecourseModel1.solve(heuristicSolution);
    }

    // Run the high-dimensional model.
    SVRPSolution solution;
    CutPool *cutPool = NULL;
    if (params.sriCuts) {
        cutPool = new CutPool();
        newParams.timeLimit = 60;
        newParams.basicVRPSD = params.basicVRPSD;
        NodeScenariosModel nodeScenariosModel(instance, newParams, *cutPool);
        if (heuristicSolution.upperBound < 1e6) {
            nodeScenariosModel.setPrimalSolution(heuristicSolution);
        }
        nodeScenariosModel.solve(solution);
        cutPool->rootBound = solution.rootBound;
    }

    // Solve the projected model.
    NodeRecourseModel nodeRecourseModel2(instance, params, cutPool);
    if (heuristicSolution.upperBound < 1e6) {
        nodeRecourseModel2.setWarmStart(heuristicSolution);
    }
    nodeRecourseModel2.solve(solution);

    if (cutPool != NULL) {
        delete cutPool;
    }

    auto done = chrono::high_resolution_clock::now();
    solution.time =
        static_cast<double>(
            chrono::duration_cast<chrono::milliseconds>(done - started)
                .count()) /
        1000.0;

    solution.print();

    // Create output as desired.
    if (params.showGraph) {
        graphviewer::viewVRPSolution(instance, solution, "SVRP solution");
    } else {
        if (params.outputFile.size() > 0) {
            solution.save(params.outputFile);
        }
        if (params.tableName.size() > 0) {
            filehandler::addEntryToTable(params, instance, solution);
        }
    }

    return 0;
}

Params getParams(int argc, const char *argv[]) {
    Params params;

    try {
        cxxopts::Options options("VRPSD", "VRPSD Program Options");
        options.add_options()("help", "produce help message")(
            "i,input", "Input file (required)",
            cxxopts::value<std::string>(params.inputFile))(
            "o,output", "Output file",
            cxxopts::value<std::string>(params.outputFile))(
            "table", "Table file",
            cxxopts::value<std::string>(params.tableName))(
            "g,graph", "Create visualization",
            cxxopts::value<bool>(params.showGraph)->default_value("false"))(
            "t,timelimit", "Time limit",
            cxxopts::value<int>(params.timeLimit)->default_value("1800"))(
            "k,number_routes", "Number of routes",
            cxxopts::value<int>(params.k))(
            "readk", "Read number of routes from input file",
            cxxopts::value<bool>(params.readKfromFile)->default_value("true"))(
            "sri", "Use SRI cuts.",
            cxxopts::value<bool>(params.sriCuts)->default_value("false"))(
            "pr", "Use partial route cuts.",
            cxxopts::value<bool>(params.partialRouteCuts)
                ->default_value("false"))(
            "ps", "Separate Parada set cuts.",
            cxxopts::value<bool>(params.paradaSet)->default_value("false"))(
            "scenarioOptimal", "Use scenario optimal recourse policy.",
            cxxopts::value<bool>()->default_value("false"))(
            "scenarioOptimalPRCuts",
            "Separate scenario optimal partial route cuts.",
            cxxopts::value<bool>(params.scenarioOptimalPRCuts)
                ->default_value("false"))(
            "basic",
            "Deactivate fixed number of routes and expected capacity "
            "constraints.",
            cxxopts::value<bool>(params.basicVRPSD)->default_value("false"))(
            "flow", "Use max-flow heuristic to separate aggregated SRIs.",
            cxxopts::value<bool>(params.sriFlowSeparation)
                ->default_value("false"))(
            "inout", "Use in and out separation.",
            cxxopts::value<bool>(params.inOut)->default_value("false"));

        auto result = options.parse(argc, argv);

        if (result.count("help") || result.count("input") == 0) {
            std::cout << options.help() << std::endl;
            exit(0);
        }

        if (params.k < 1 && !params.readKfromFile) {
            std::cout << "Invalid value of k and not set to read k from file."
                      << std::endl;
            exit(1);
        }

        if (result["scenarioOptimal"].as<bool>()) {
            params.policy = SCENARIO_OPTIMAL;
            params.scenarioOptimalPRCuts = true;
        }
    } catch (std::exception &e) {
        std::cerr << "Program Options Error: " << e.what() << "\n";
        exit(1);
    } catch (...) {
        std::cerr << "Unknown error!" << "\n";
        exit(1);
    }

    if (params.scenarioOptimalPRCuts && !params.partialRouteCuts) {
        std::cerr << "Error: "
                  << "improved partial routes can only be used when partial "
                     "routes are used."
                  << std::endl;
        exit(1);
    }

    if (params.sriFlowSeparation && !params.sriCuts) {
        std::cerr << "Error: max-flow heuristic can only be used when "
                     "separating aggregated SRIs."
                  << std::endl;
        exit(1);
    }

    return params;
}