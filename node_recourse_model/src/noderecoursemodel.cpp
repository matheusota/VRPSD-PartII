#include "noderecoursemodel.h"

NodeRecourseModel::NodeRecourseModel(SVRPInstance &instance,
                                     const Params &params, CutPool *cutPool_)
    : instance(instance), params(params), model(env), x(instance.g),
      y(instance.g), cutBuilder(instance, params, x, y),
      cvrpsepSeparator(instance, params), paradaSeparator(instance, params),
      partialRouteSeparator(instance, params),
      aggregatedSRISeparator(instance, params),
      callback(instance, params, x, y, cutBuilder, cvrpsepSeparator,
               paradaSeparator, partialRouteSeparator, aggregatedSRISeparator) {
    cutPool = cutPool_;

    // Set model.
    setBasicModel();

    // Set Gurobi callback.
    model.setCallback(&callback);
}

void NodeRecourseModel::setWarmStart(const SVRPSolution &primalSolution) {
    std::vector<Edge> solutionEdges =
        primalSolution.getEdgesInSolution(instance);
    EdgeValueMap xValue(instance.g, 0.0);
    NodeValueMap yValue(instance.g, 0.0);
    ClassicalRecourseHelper classicalRecourseHelper(instance);

    for (Edge e : solutionEdges) {
        xValue[e] += 1.0;
    }

    for (const auto &route : primalSolution.routes) {
        NodeVectorValueMap failures(
            instance.g, std::vector<double>(instance.nScenarios, 0.0));
        classicalRecourseHelper.getRouteRecourseCostAndFailures(route,
                                                                failures);

        for (int id : route) {
            Node v = instance.g.nodeFromId(id);

            for (int scenarioId = 0; scenarioId < instance.nScenarios;
                 scenarioId++) {
                yValue[v] += (instance.getEdgeRecourseCost(v) *
                              failures[v][scenarioId]) /
                             static_cast<double>(instance.nScenarios);
            }
        }
    }

    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        x[e].set(GRB_DoubleAttr_Start, xValue[e]);
    }

    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) != instance.depot) {
            y[v].set(GRB_DoubleAttr_Start, yValue[v]);
        }
    }
}

void NodeRecourseModel::setBasicModel() {
    model.set(GRB_StringAttr_ModelName, "GurobiModel");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Must set to use lazy contraints.
    model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
    model.getEnv().set(GRB_IntParam_PreCrush, 1);

    // Some other parameters.
    model.getEnv().set(GRB_DoubleParam_TimeLimit, params.timeLimit);
    model.getEnv().set(GRB_IntParam_Threads, 1);

    // Set x variables.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        double UB = 1.0;
        if (instance.g.id(instance.g.u(e)) == instance.depot ||
            instance.g.id(instance.g.v(e)) == instance.depot) {
            UB = 2.0;
        }

        x[e] = model.addVar(
            0.0, UB, instance.weight[e], GRB_INTEGER,
            "x_" + std::to_string(instance.g.id(instance.g.u(e))) + "," +
                std::to_string(instance.g.id(instance.g.v(e))));
    }

    // Set y variables.
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) != instance.depot) {
            y[v] = model.addVar(0.0, instance.nScenarios,
                                instance.getEdgeRecourseCost(v) /
                                    static_cast<double>(instance.nScenarios),
                                GRB_CONTINUOUS,
                                "y_" + std::to_string(instance.g.id(v)));
        }
    }
    model.update();

    // Incidence constraints.
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        int id = instance.g.id(v);
        double RHS = (id == instance.depot) ? 2.0 * instance.k : 2.0;

        GRBLinExpr expr = 0.0;
        for (IncEdgeIt e(instance.g, v); e != INVALID; ++e) {
            expr += x[e];
        }

        if (id == instance.depot && params.basicVRPSD) {
            model.addConstr(expr >= 2.0, "incidence_" + std::to_string(id));
        } else {
            model.addConstr(expr == RHS, "incidence_" + std::to_string(id));
        }
    }
    model.update();
}

bool NodeRecourseModel::solve(SVRPSolution &solution) {
    double MIPTimeLimit = params.timeLimit;
    auto started = chrono::high_resolution_clock::now();

    // Set initial root bound to check later when using Lagrangian cuts.
    initialRootBound = solution.rootBound;

    // Solve the root separately.
    // First we change the variable types.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        x[e].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
    }
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) != 0) {
            y[v].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
        }
    }
    model.update();

    // Warm start with some cuts.
    if (cutPool != nullptr) {
        const std::vector<CutData> initialCuts = cutPool->getAllCuts();
        addInitialCuts(initialCuts);
    }

    // Solve root and get LP bound.
    solveRootLP(started);
    solution.rootBound = model.get(GRB_DoubleAttr_ObjVal);

    // Change variable types back to integer.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        x[e].set(GRB_CharAttr_VType, GRB_INTEGER);
    }

    // Delete inactive cuts.
    // cleanupCuts();

    auto done = chrono::high_resolution_clock::now();
    solution.rootTime +=
        static_cast<double>(
            chrono::duration_cast<chrono::milliseconds>(done - started)
                .count()) /
        1000.0;
    MIPTimeLimit -= solution.rootTime;

    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Root bound: " << solution.rootBound << std::endl;
    std::cout << "Root time: " << solution.rootTime << " s" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;

    if (MIPTimeLimit < 0) {
        return 0;
    }

    // Solve the MIP.
    model.getEnv().set(GRB_DoubleParam_TimeLimit, MIPTimeLimit);
    cvrpsepSeparator.setMaxNoOfCapCuts(2);
    model.update();
    try {
        model.optimize();
    } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Error during optimization" << endl;
    }

    // Collect solution and statistics.
    setSolution(solution);

    return 0;
}

void NodeRecourseModel::setSolution(SVRPSolution &solution) {
    solution.lowerBound = model.get(GRB_DoubleAttr_ObjBound);
    solution.upperBound = model.get(GRB_DoubleAttr_ObjVal);
    solution.nodesExplored = model.get(GRB_DoubleAttr_NodeCount);
    solution.gap = model.get(GRB_DoubleAttr_MIPGap);
    solution.rootBound = std::max(solution.rootBound, callback.rootBound);
    solution.cvrpsepCuts += cvrpsepSeparator.nCuts;
    solution.cvrpsepTime += cvrpsepSeparator.time;
    solution.paradaSetCuts += paradaSeparator.nSetCuts;
    solution.paradaTime += paradaSeparator.time;
    solution.partialRouteCuts += partialRouteSeparator.nCuts;
    solution.partialRouteTime += partialRouteSeparator.time;
    solution.aggregatedSriCuts += aggregatedSRISeparator.nCuts;
    solution.aggregatedSriTime += aggregatedSRISeparator.time;

    // Get back the solution.
    if (solution.upperBound <= 1e6) {
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            solution.solved = true;
        }

        solution.deterministicCost = 0.0;
        solution.modelRecourseCost = 0.0;

        EdgeIntMap edgeCount(instance.g, 0);
        for (EdgeIt e(instance.g); e != INVALID; ++e) {
            if (x[e].get(GRB_DoubleAttr_X) > 1e-4) {
                edgeCount[e] =
                    static_cast<int>(std::round(x[e].get(GRB_DoubleAttr_X)));
                std::cout << x[e].get(GRB_StringAttr_VarName) << " = "
                          << x[e].get(GRB_DoubleAttr_X) << std::endl;
                solution.deterministicCost += edgeCount[e] * instance.weight[e];
            }
        }
        solution.setRoutesFromSolution(instance, edgeCount);

        for (NodeIt v(instance.g); v != INVALID; ++v) {
            if (instance.g.id(v) != instance.depot &&
                y[v].get(GRB_DoubleAttr_X) > 1e-4) {
                std::cout << y[v].get(GRB_StringAttr_VarName) << " = "
                          << y[v].get(GRB_DoubleAttr_X) << std::endl;
                solution.modelRecourseCost +=
                    (y[v].get(GRB_DoubleAttr_X) *
                     instance.getEdgeRecourseCost(v)) /
                    static_cast<double>(instance.nScenarios);
            }
        }

        solution.cost = solution.deterministicCost + solution.modelRecourseCost;
        solution.recourseCost =
            instance.classicalRecourseHelper.getRecourseCost(solution.routes);
        solution.optimalRecourseCost =
            instance.optimalRecourseHelper.getRecourseCost(solution.routes);

        std::cout << "Deterministic cost " << solution.deterministicCost
                  << endl;
        std::cout << "Model recourse cost " << solution.modelRecourseCost
                  << endl;
        std::cout << "Optimal recourse cost " << solution.optimalRecourseCost
                  << endl;
        std::cout << "Classical recourse cost " << solution.recourseCost
                  << endl;
        std::cout << "Root bound: " << solution.rootBound << std::endl;

        // We check the recourse cost only when the problem was solved.
        if (solution.solved) {
            if (params.policy == SCENARIO_OPTIMAL) {
                assert(std::abs(solution.modelRecourseCost -
                                solution.optimalRecourseCost) <= 1e-4);
            } else {
                assert(std::abs(solution.modelRecourseCost -
                                solution.recourseCost) <= 1e-4);
            }
        }
    } else {
        std::cout << "Did not find optimal solution." << std::endl;
    }
}

void NodeRecourseModel::solveRootLP(
    std::chrono::time_point<std::chrono::high_resolution_clock> &started) {
    double prevDualBound = 0.0;
    int noImprovCount = 0;
    bool solveLP = true;
    EdgeValueMap xValue(instance.g, 0.0);
    NodeValueMap yValue(instance.g, 0.0);
    std::vector<CutData> separatedCuts;
    bool addedLagrangianCut = false;

    while (solveLP && noImprovCount <= 10) {
        solveLP = false;
        model.update();
        model.optimize();

        auto done = chrono::high_resolution_clock::now();
        double elapsedTime =
            static_cast<double>(
                chrono::duration_cast<chrono::milliseconds>(done - started)
                    .count()) /
            1000.0;
        if (elapsedTime > params.timeLimit) {
            return;
        }

        // Get dual bound and update no improvement counter.
        double currDualBound = model.get(GRB_DoubleAttr_ObjBound);
        if (currDualBound - prevDualBound >= 1e-6) {
            noImprovCount = 0;
        } else {
            noImprovCount++;
        }
        prevDualBound = currDualBound;

        // Construct current solution.
        double totalRecourse;
        getCurrentSolution(xValue, yValue);

        // Separate additional cuts.
        callback.separateCuts(xValue, yValue, separatedCuts, false, 0);

        // Add Lagrangian cut if necessary.
        if (separatedCuts.empty() && cutPool != nullptr &&
            currDualBound <= cutPool->rootBound - 0.01) {
            assert(!addedLagrangianCut);
            const std::vector<CutData> &cutPoolCuts = cutPool->getAllCuts();
            assert(cutPoolCuts.back().name == "LagrangianCut");
            separatedCuts.push_back(cutPoolCuts.back());
            addedLagrangianCut = true;
        }

        if (addSeparatedCuts(separatedCuts)) {
            solveLP = true;
        }
    }

    separatedCuts.clear();
}

void NodeRecourseModel::getCurrentSolution(EdgeValueMap &xValue,
                                           NodeValueMap &yValue) {
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        xValue[e] = x[e].get(GRB_DoubleAttr_X);
    }

    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) != instance.depot) {
            yValue[v] = y[v].get(GRB_DoubleAttr_X);
        }
    }
}

bool NodeRecourseModel::addSeparatedCuts(
    const std::vector<CutData> &separatedCuts) {
    bool addedCut = false;

    for (const CutData &cutData : separatedCuts) {
        assert(std::isfinite(cutData.LHS));
        assert(std::isfinite(cutData.RHS));
        double sign = (cutData.sense == GRB_LESS_EQUAL) ? 1 : -1;
        if (sign * (cutData.LHS - cutData.RHS) <= 1e-6) {
            continue;
        }

        try {
            GRBConstr addedConstraint =
                model.addConstr(cutBuilder.buildCutExpr(cutData), cutData.sense,
                                cutData.RHS, cutData.name);
            addedCut = true;
        } catch (GRBException e) {
            std::cout << "Gurobi exception: " << e.getMessage()
                      << " code: " << e.getErrorCode() << std::endl;
            exit(1);
        }
    }

    return addedCut;
}

void NodeRecourseModel::addInitialCuts(
    const std::vector<CutData> &initialCuts) {
    std::unordered_set<std::vector<int>, hash_vector> addedSRI;
    std::unordered_set<std::vector<int>, hash_vector> addedCVRPSEP;

    for (const CutData &cutData : initialCuts) {
        if (cutData.name == "CVRPSEP") {
            if (addedCVRPSEP.find(cutData.customers) != addedCVRPSEP.end()) {
                continue;
            } else {
                addedCVRPSEP.insert(cutData.customers);
                GRBConstr addedConstraint =
                    model.addConstr(cutBuilder.buildCutExpr(cutData),
                                    cutData.sense, cutData.RHS, cutData.name);
            }
        } else if (cutData.name == "SRI" && std::abs(cutData.dual) >= 1e-6) {
            if (addedSRI.find(cutData.customers) != addedSRI.end()) {
                continue;
            } else {
                addedSRI.insert(cutData.customers);

                // // We collect these scenarios in the aggregated sri
                // separator. std::vector<int> scenarioIds; for (auto &p :
                // cutData.scenarioPairs) {
                //     scenarioIds.push_back(p.first);
                // }
                // aggregatedSRISeparator.addCriticalScenarios(scenarioIds);

                // Do some ugly stuff to add the aggregated set cut.
                std::vector<CutData> tmpCuts;
                EdgeValueMap tmpEdgeMap(instance.g, 1.0);
                NodeValueMap tmpNodeMap(instance.g, 0.0);
                aggregatedSRISeparator.addCutFromSet(
                    tmpEdgeMap, tmpNodeMap, cutData.customers, tmpCuts);

                for (CutData &tmpData : tmpCuts) {
                    GRBConstr addedConstraint = model.addConstr(
                        cutBuilder.buildCutExpr(tmpData), tmpData.sense,
                        tmpData.RHS, tmpData.name);
                    // addedConstraint.set(GRB_IntAttr_Lazy, 1);
                }
            }
        }
    }
}

void NodeRecourseModel::cleanupCuts() {
    GRBConstr *constraints = model.getConstrs();
    int numConstrs = model.get(GRB_IntAttr_NumConstrs);

    std::vector<GRBConstr> toRemove;
    int removedCuts = 0;

    for (int i = 0; i < numConstrs; ++i) {
        std::string name = constraints[i].get(GRB_StringAttr_ConstrName);

        if (name.find("incidence") == std::string::npos) {

            double slack = constraints[i].get(GRB_DoubleAttr_Slack);

            // If slack is large, the cut is not active
            if (std::abs(slack) >= 0.9) {
                toRemove.push_back(constraints[i]);
            }
        }
    }

    // Remove the inactive cuts
    for (auto &constr : toRemove) {
        model.remove(constr);
        removedCuts++;
    }

    delete[] constraints;
    model.update();

    std::cout << "Removed " << removedCuts << " inactive cuts." << std::endl;
}