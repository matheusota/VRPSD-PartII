#include "filehandler.h"

namespace filehandler {
bool readInstance(const Params &params, SVRPInstance &instance) {
    std::ifstream file;
    std::string line;
    std::vector<std::string> lines;
    std::vector<double> demands;
    bool hasDistanceMatrixSection = false;

    file.open(params.inputFile.c_str());
    if (!file) {
        cerr << "File '" << params.inputFile << "' does not exist.\n";
        exit(0);
    }

    // put file lines in an array
    while (std::getline(file, line)) {
        lines.push_back(line);
    }

    // read each line
    size_t i = 0;
    while (i < lines.size()) {
        line = lines[i];

        if (line.find("NAME") != std::string::npos) {
            std::smatch match;
            regex_search(line, match, regex("n([0-9]+)-k([0-9]+)"));

            if (match.str(2) != "" && params.readKfromFile) {
                instance.k = stoi(match.str(2));
            }
        } else if (line.find("DIMENSION") != std::string::npos) {
            std::smatch match;
            regex_search(line, match, regex("[0-9]+"));
            instance.n = stoi(match.str(0));
            instance.m = instance.n * (instance.n - 1);
            demands.resize(instance.n);

            if (!params.readKfromFile) {
                instance.k = params.k;
            }
        } else if (line.find("NUMBER_OF_SCENARIOS") != std::string::npos) {
            std::smatch match;
            regex_search(line, match, regex("[0-9]+"));
            instance.nScenarios = stoi(match.str(0));
        } else if (line.find("CAPACITY") != std::string::npos) {
            std::smatch match;
            regex_search(line, match, regex("[0-9]+"));
            instance.capacity = stoi(match.str(0));
        } else if (line.find("N_VEHICLES") != std::string::npos) {
            if (params.readKfromFile) {
                std::smatch match;
                regex_search(line, match, regex("[0-9]+"));
                instance.k = stoi(match.str(0));
            }
        } else if (line.find("NODE_COORD_SECTION") != std::string::npos) {
            std::smatch match;
            line = lines[++i];
            while (std::regex_search(line, match,
                                     regex("\\s*([-+]?[0-9]+)\\s*([-+]?[0-9]+)"
                                           "\\s*([-+]?[0-9]+)\\s*"))) {
                int id = stoi(match.str(1));
                double posx = stof(match.str(2));
                double posy = stof(match.str(3));

                Node v = instance.g.addNode();
                instance.posx[v] = posx;
                instance.posy[v] = posy;

                DNode d_v = instance.d_g.addNode();
                instance.d_posx[d_v] = posx;
                instance.d_posy[d_v] = posy;

                line = lines[++i];
            }
            i--;
        }
        // We use this section only if we are not using scenarios.
        else if (line.find("DEMAND_SECTION") != std::string::npos) {
            std::smatch match;
            int id;
            double demand;

            line = lines[++i];
            while (std::regex_search(
                line, match, regex("\\s*([0-9]+)\\s*([0-9]+\\.?[0-9]*)\\s*"))) {
                id = stoi(match.str(1));

                demand = stof(match.str(2));

                Node v = instance.g.nodeFromId(id - 1);
                instance.demand[v] = demand;

                DNode d_v = instance.d_g.nodeFromId(id - 1);
                instance.d_demand[d_v] = demand;

                line = lines[++i];
            }
            i--;
        } else if (line.find("DISTANCE_MATRIX") != std::string::npos) {
            hasDistanceMatrixSection = true;
            line = lines[++i];
            double roundFactor = 1.0;

            for (int j = 0; j < instance.n; j++) {
                stringstream ss(line);
                string s;
                int k = 0;

                Node source = instance.g.nodeFromId(j);
                DNode d_source = instance.d_g.nodeFromId(j);
                while (std::getline(ss, s, ' ')) {
                    Node target = instance.g.nodeFromId(k);
                    DNode d_target = instance.d_g.nodeFromId(k);
                    if (k < instance.n) {
                        if (j < k) {
                            Edge e = instance.g.addEdge(source, target);
                            instance.weight[e] =
                                std::round(stod(s) * roundFactor) / roundFactor;
                        }
                        if (j != k) {
                            Arc a = instance.d_g.addArc(d_source, d_target);
                            instance.d_weight[a] =
                                std::round(stod(s) * roundFactor) / roundFactor;
                        }
                    }
                    k++;
                }

                line = lines[++i];
            }
            i--;
        } else if (line.find("DEMAND_SCENARIO_SECTION") != std::string::npos) {
            // allocate scenarios matrix
            instance.scenariosMatrix.resize(
                instance.n + 1, std::vector<int>(instance.nScenarios, 0));

            // fill scenarios matrix
            line = lines[++i];
            for (int j = 0; j < instance.n; j++) {
                stringstream ss(line);
                string s;
                int k = 0;

                Node v = instance.g.nodeFromId(j);
                DNode d_v = instance.d_g.nodeFromId(j);
                double accDemand = 0.0;
                while (std::getline(ss, s, ' ')) {
                    if (k < instance.nScenarios) {
                        instance.scenariosMatrix[j][k] = stoi(s);
                        assert(instance.scenariosMatrix[j][k] <=
                               static_cast<int>(instance.capacity));
                        accDemand += instance.scenariosMatrix[j][k];
                    }
                    k++;
                }

                // set the demand as the average
                instance.demand[v] =
                    std::round(accDemand / instance.nScenarios);
                instance.d_demand[d_v] =
                    std::round(accDemand / instance.nScenarios);
                line = lines[++i];
                demands[j] = instance.demand[v];
            }
            i--;
        }

        // go to next line
        i++;
    }

    // add edges (it is always a complete graph)
    instance.m = (instance.n * (instance.n - 1)) / 2;
    if (!hasDistanceMatrixSection) {
        for (NodeIt v(instance.g); v != INVALID; ++v) {
            for (NodeIt u(instance.g); u != INVALID; ++u) {
                if (instance.g.id(v) < instance.g.id(u)) {
                    Edge e = instance.g.addEdge(v, u);
                    instance.weight[e] = std::round(
                        hypot(abs(instance.posx[v] - instance.posx[u]),
                              abs(instance.posy[v] - instance.posy[u])));
                }
            }
        }

        for (DNodeIt v(instance.d_g); v != INVALID; ++v) {
            for (DNodeIt u(instance.d_g); u != INVALID; ++u) {
                if (instance.d_g.id(v) != instance.d_g.id(u)) {
                    Arc a = instance.d_g.addArc(v, u);
                    instance.d_weight[a] = std::round(
                        hypot(abs(instance.d_posx[v] - instance.d_posx[u]),
                              abs(instance.d_posy[v] - instance.d_posy[u])));
                }
            }
        }
    }

    // Sort scenarios by their total demands.
    instance.sortScenarios();

    // Set upper bound map.
    double totalDemand = 0.0;
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        totalDemand += instance.demand[v];
    }

    for (NodeIt v(instance.g); v != INVALID; ++v) {
        double tmpDemand = totalDemand - instance.demand[v];
        instance.upperBoundMap[v] =
            std::ceil(tmpDemand / instance.capacity) <= instance.k - 1.0 ? 2
                                                                         : 1;
    }

    // depot
    instance.depot = 0;

    return true;
} // namespace filehandler

bool saveSolution(const SVRPSolution &sol, std::string outputName) {
    std::ofstream myfile;
    myfile.open(outputName);

    if (!myfile) {
        cerr << "Error writing to " << outputName << endl;
        exit(0);
    }

    for (auto route : sol.routes) {
        for (int i : route) {
            myfile << to_string(i) << " ";
        }
        myfile << endl;
    }

    myfile.close();
    return true;
}

bool addEntryToTable(const Params &params, const SVRPInstance &instance,
                     const SVRPSolution &solution) {
    std::ofstream data;
    data.open(params.tableName, std::ios_base::app);

    if (!data) {
        cerr << "Error writing to " << params.tableName << endl;
        exit(0);
    }

    std::string algorithm = "SVRP";
    if (params.gendreauCuts) {
        algorithm += "-Gendreau";
    }
    if (params.sriCuts) {
        algorithm += "-SRI";
    }
    if (params.paradaSet) {
        algorithm += "-PARADA_SET";
    }
    if (params.partialRouteCuts) {
        algorithm += "-PARTIAL_ROUTE";
    }
    if (params.alg == NODE) {
        algorithm += "-NODE";
    }
    if (params.alg == NODE_SCENARIOS) {
        algorithm += "-NODE+SCEN";
    }
    if (params.mipSeparation) {
        algorithm += "-MIPHEUR";
    }
    if (params.benders) {
        algorithm += "-BENDERS";
    }
    if (params.lagrangian) {
        algorithm += "-LAG";
    }
    if (params.scenarioOptimalPRCuts) {
        algorithm += "-PR_SCENOPT_CUTS";
    }
    if (params.basicVRPSD) {
        algorithm += "-Basic";
    }
    if (params.inOut) {
        algorithm += "-InOut";
    }
    if (params.sriFlowSeparation) {
        algorithm += "-FLOWSEP";
    }
    if (params.dualSetCut) {
        algorithm += "-DUALSET";
    }

    data << algorithm << " , " << params.inputFile << " , " << solution.cost
         << " , " << solution.rootBound << " , "
         << (std::abs(solution.cost - solution.rootBound) / solution.cost) *
                100.0
         << " , " << solution.rootTime << " , " << solution.time << " , "
         << instance.k << " , " << solution.deterministicCost << " , "
         << solution.modelRecourseCost << " , " << solution.optimalRecourseCost
         << " , " << solution.recourseCost << " , " << solution.upperBound
         << " , " << solution.lowerBound << " , " << solution.gap * 100.0
         << " , " << solution.nodesExplored << " , " << solution.cvrpsepCuts
         << " , " << solution.cvrpsepTime << " , " << solution.partialRouteCuts
         << " , " << solution.partialRouteTime << " , "
         << solution.paradaPathCuts << " , " << solution.paradaSetCuts << " , "
         << solution.paradaTime << " , " << solution.sriCuts << " , "
         << solution.sriTime << " , " << solution.aggregatedSriCuts << " , "
         << solution.aggregatedSriTime << " , " << (solution.solved ? 1 : 0)
         << std::endl;

    data.close();
    return true;
}
} // namespace filehandler
