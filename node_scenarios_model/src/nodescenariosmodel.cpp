#include "nodescenariosmodel.h"

NodeScenariosModel::NodeScenariosModel(const SVRPInstance &instance,
                                       const Params &params, CutPool *cutPool_)
    : instance(instance), params(params), model(env), x(instance.g),
      y(instance.g), interiorX(instance.g, 0.0), interiorY(instance.g),
      cvrpsepSeparator(instance, params),
      sriSeparator(instance, params, cvrpsepSeparator),
      callback(instance, params, x, y, cvrpsepSeparator, sriSeparator) {
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        interiorY[v].resize(instance.nScenarios);
    }

    // Set model.
    setBasicModel();
    cutPool = cutPool_;
    model.setCallback(&callback);
}

void NodeScenariosModel::setPrimalSolution(const SVRPSolution &primalSolution) {
    UB = 0.0;
    std::vector<Edge> solutionEdges =
        primalSolution.getEdgesInSolution(instance);
    ClassicalRecourseHelper classicalRecourseHelper(instance);

    for (Edge e : solutionEdges) {
        interiorX[e] += 1;
        UB += instance.weight[e];
    }

    for (auto &route : primalSolution.routes) {
        NodeVectorValueMap failures(
            instance.g, std::vector<double>(instance.nScenarios, 0.0));
        classicalRecourseHelper.getRouteRecourseCostAndFailures(route,
                                                                failures);

        for (int id : route) {
            Node v = instance.g.nodeFromId(id);

            for (int scenarioId = 0; scenarioId < instance.nScenarios;
                 scenarioId++) {
                interiorY[v][scenarioId] = failures[v][scenarioId];
                UB += (instance.getEdgeRecourseCost(v) *
                       failures[v][scenarioId]) /
                      static_cast<double>(instance.nScenarios);
            }
        }
    }
}

void NodeScenariosModel::setBasicModel() {
    model.set(GRB_StringAttr_ModelName, "GurobiModel");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Must set to use lazy contraints.
    model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
    model.getEnv().set(GRB_IntParam_PreCrush, 1);

    // Some other parameters.
    model.getEnv().set(GRB_DoubleParam_TimeLimit, params.timeLimit);
    model.getEnv().set(GRB_IntParam_Threads, 1);
    model.getEnv().set(GRB_DoubleParam_MIPGap, 1e-9);
    model.getEnv().set(GRB_IntParam_Cuts, 0);
    model.getEnv().set(GRB_IntParam_Presolve, 2);
    model.getEnv().set(GRB_IntParam_Method, 1);

    // Set x variables.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        x[e] = model.addVar(
            0.0, GRB_INFINITY, instance.weight[e], GRB_INTEGER,
            "x_" + std::to_string(instance.g.id(instance.g.u(e))) + "," +
                std::to_string(instance.g.id(instance.g.v(e))));
    }

    // Set y variables.
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) != instance.depot) {
            for (int scenarioId = 0; scenarioId < instance.nScenarios;
                 scenarioId++) {
                y[v].push_back(model.addVar(
                    0.0, GRB_INFINITY,
                    (1.0 / static_cast<double>(instance.nScenarios)) *
                        instance.getEdgeRecourseCost(v),
                    GRB_CONTINUOUS,
                    "y_" + std::to_string(instance.g.id(v)) + "^" +
                        std::to_string(scenarioId)));
            }
        }
    }

    // Bound constraints.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        double UB = 1.0;
        if (instance.g.id(instance.g.u(e)) == instance.depot ||
            instance.g.id(instance.g.v(e)) == instance.depot) {
            UB = 2.0;
        }

        model.addConstr(x[e] <= UB,
                        "edgebound_" + std::to_string(instance.g.id(e)) + "_" +
                            std::to_string(instance.g.id(instance.g.u(e))) +
                            "," +
                            std::to_string(instance.g.id(instance.g.v(e))));
    }

    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) != instance.depot) {
            for (int scenarioId = 0; scenarioId < instance.nScenarios;
                 scenarioId++) {
                model.addConstr(y[v][scenarioId] <= 1.0,
                                "nodebound_" +
                                    std::to_string(instance.g.id(v)) + "_" +
                                    std::to_string(scenarioId));
            }
        }
    }

    // Incidence constraints.
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        int v_id = instance.g.id(v);
        double RHS = (v_id == instance.depot) ? 2.0 * instance.k : 2.0;

        GRBLinExpr expr = 0.0;
        for (IncEdgeIt e(instance.g, v); e != INVALID; ++e) {
            expr += x[e];
        }

        if (v_id == instance.depot && params.basicVRPSD) {
            model.addConstr(expr >= 2.0, "incidence_" + std::to_string(v_id));
        } else {
            model.addConstr(expr == RHS, "incidence_" + std::to_string(v_id));
        }
    }
    model.update();
}

bool NodeScenariosModel::solve(SVRPSolution &solution) {
    double MIPTimeLimit = params.timeLimit;
    auto started = chrono::high_resolution_clock::now();

    // Solve the root separately.
    // Change variable types.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        x[e].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
    }
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) != instance.depot) {
            for (int scenarioId = 0; scenarioId < instance.nScenarios;
                 scenarioId++) {
                y[v][scenarioId].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
            }
        }
    }
    model.update();

    // Solve root.
    cvrpsepSeparator.setMaxNoOfCapCuts(10);
    solveRootLP(started);
    solution.rootBound = model.get(GRB_DoubleAttr_ObjVal);

    // Change variable types back to integer.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        x[e].set(GRB_CharAttr_VType, GRB_INTEGER);
    }

    auto done = chrono::high_resolution_clock::now();
    solution.rootTime =
        static_cast<double>(
            chrono::duration_cast<chrono::milliseconds>(done - started)
                .count()) /
        1000.0;
    MIPTimeLimit -= solution.rootTime;

    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Root bound: " << solution.rootBound << std::endl;
    std::cout << "Root time: " << solution.rootTime << " s" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;

    // If using cut pool, build dual cut and return.
    if (cutPool != nullptr) {
        solution.cvrpsepCuts += cvrpsepSeparator.nCuts;
        solution.cvrpsepTime += cvrpsepSeparator.time;
        solution.sriCuts += sriSeparator.nCuts;
        solution.sriTime += sriSeparator.time;

        double dualObj = buildLagrangianCutAndGetDualObj();
        assert(std::abs(dualObj - solution.rootBound) <= 1e-4);
        return 0;
    }

    if (MIPTimeLimit < 0) {
        return 0;
    }

    // Solve the MIP.
    model.getEnv().set(GRB_DoubleParam_TimeLimit, MIPTimeLimit);
    // cvrpsepSeparator->setMaxNoOfCapCuts(5);
    model.update();
    model.optimize();

    // Collect solution and statistics.
    setSolution(solution);

    return 0;
}

void NodeScenariosModel::setSolution(SVRPSolution &solution) {
    solution.lowerBound = model.get(GRB_DoubleAttr_ObjBound);
    solution.upperBound = model.get(GRB_DoubleAttr_ObjVal);
    solution.nodesExplored = model.get(GRB_DoubleAttr_NodeCount);
    solution.gap = model.get(GRB_DoubleAttr_MIPGap);
    // solution.rootBound = std::max(solution.rootBound, callback->rootBound);
    solution.cvrpsepCuts += cvrpsepSeparator.nCuts;
    solution.cvrpsepTime += cvrpsepSeparator.time;
    solution.sriCuts += sriSeparator.nCuts;
    solution.sriTime += sriSeparator.time;
    solution.deterministicCost = 0.0;
    solution.recourseCost = 0.0;
    solution.modelRecourseCost = 0.0;

    // Get back the solution.
    if (solution.upperBound <= 1e6) {
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            solution.solved = true;
        }

        EdgeIntMap edgeCount(instance.g, 0);
        for (EdgeIt e(instance.g); e != INVALID; ++e) {
            if (x[e].get(GRB_DoubleAttr_X) > 0.001) {
                edgeCount[e] =
                    static_cast<int>(std::round(x[e].get(GRB_DoubleAttr_X)));
                std::cout << x[e].get(GRB_StringAttr_VarName) << " = "
                          << x[e].get(GRB_DoubleAttr_X) << std::endl;
                solution.deterministicCost += edgeCount[e] * instance.weight[e];
            }
        }
        solution.setRoutesFromSolution(instance, edgeCount);

        for (int scenarioId = 0; scenarioId < instance.nScenarios;
             scenarioId++) {
            for (NodeIt v(instance.g); v != INVALID; ++v) {
                if (instance.g.id(v) != instance.depot &&
                    y[v][scenarioId].get(GRB_DoubleAttr_X) > 0.001) {
                    std::cout
                        << y[v][scenarioId].get(GRB_StringAttr_VarName) << " = "
                        << y[v][scenarioId].get(GRB_DoubleAttr_X)
                        << " recourse cost: " << instance.getEdgeRecourseCost(v)
                        << std::endl;
                    solution.modelRecourseCost +=
                        (1.0 / static_cast<double>(instance.nScenarios)) *
                        instance.getEdgeRecourseCost(v) *
                        y[v][scenarioId].get(GRB_DoubleAttr_X);
                }
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
    } else {
        std::cout << "Didn’t find optimal solution." << std::endl;
    }
}

// This solves the LP and update the out point and the candidate solution.
void NodeScenariosModel::solveLPAndUpdateSolutions(EdgeValueMap &xValue,
                                                   NodeVectorValueMap &yValue,
                                                   EdgeValueMap &candX,
                                                   NodeVectorValueMap &candY) {
    model.update();
    model.optimize();
    assert(model.get(GRB_IntAttr_Status) == GRB_OPTIMAL);
    LB = model.get(GRB_DoubleAttr_ObjBound);

    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        xValue[e] = x[e].get(GRB_DoubleAttr_X);
    }

    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        for (NodeIt v(instance.g); v != INVALID; ++v) {
            if (instance.g.id(v) != instance.depot) {
                yValue[v][scenarioId] = y[v][scenarioId].get(GRB_DoubleAttr_X);
            }
        }
    }

    // Update candidate.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        candX[e] = (1.0 - lambda) * interiorX[e] + lambda * xValue[e];
    }

    for (NodeIt v(instance.g); v != INVALID; ++v) {
        for (int scenarioId = 0; scenarioId < instance.nScenarios;
             scenarioId++) {
            candY[v][scenarioId] = (1.0 - lambda) * interiorY[v][scenarioId] +
                                   lambda * yValue[v][scenarioId];
        }
    }
}

void NodeScenariosModel::solveRootLP(
    std::chrono::time_point<std::chrono::high_resolution_clock> &started) {
    int it = 0;
    int countCuts = 0;
    EdgeValueMap xValue(instance.g, 0.0);
    NodeVectorValueMap yValue(instance.g,
                              std::vector<double>(instance.nScenarios, 0.0));
    EdgeValueMap candX(instance.g, 0.0);
    NodeVectorValueMap candY(instance.g,
                             std::vector<double>(instance.nScenarios, 0.0));
    std::vector<CutData> separatedCuts;
    model.set(GRB_IntParam_LogToConsole, 0);
    lambda = 1.0; // This is the convex comb. parameter for in/out.

    // Initialize x and y values.
    solveLPAndUpdateSolutions(xValue, yValue, candX, candY);

    // model.set(GRB_IntParam_NumericFocus, 3);
    while ((UB - LB) / UB > 1e-2) {
        std::cout << "NodeScenariosModel LP, Current gap: " << (UB - LB) / UB
                  << std::endl;
        std::cout << "it: " << it << " cuts added: " << countCuts
                  << " LB: " << LB << " UB: " << UB << std::endl;
        it++;

        auto done = chrono::high_resolution_clock::now();
        double elapsedTime =
            static_cast<double>(
                chrono::duration_cast<chrono::milliseconds>(done - started)
                    .count()) /
            1000.0;
        if (elapsedTime > params.timeLimit) {
            return;
        }

        // Call separation on candidate solution.
        // First call CVRPSEP.
        separatedCuts.clear();
        cvrpsepSeparator.separateCVRPSEPCuts(candX, separatedCuts);
        const std::vector<std::vector<int>> &customerSets =
            cvrpsepSeparator.getCustomerSets();
        const std::vector<int> &nVehicles = cvrpsepSeparator.getNVehicles();

        // Next, separate SRIs.
        if (separatedCuts.empty()) {
            for (size_t i = 0; i < customerSets.size(); i++) {
                sriSeparator.separateSRIFromCustomerSet(
                    candX, candY, customerSets[i], separatedCuts);
            }
        }
        sriSeparator.heuristicSeparation(candX, candY, separatedCuts);

        SRI_SEPARATION_STATUS status = SRI_SEPARATION_STATUS::NO_CUTS;
        if (separatedCuts.empty()) {
            sriSeparator.mipSeparation(candX, candY, separatedCuts, status);
            if (params.inOut) {
                lambda = 0.3;
            }
        }

        int addedCuts =
            addSeparatedCuts(xValue, yValue, candX, candY, separatedCuts);
        countCuts += addedCuts;

        // If candidate is feasible, update interior solution.
        if (addedCuts == 0) {
            // MIP separation was successful, so candidate is feasible.
            if (status == SRI_SEPARATION_STATUS::NO_CUTS) {
                UB = 0.0;

                for (EdgeIt e(instance.g); e != INVALID; ++e) {
                    interiorX[e] = candX[e];
                    UB += candX[e] * instance.weight[e];
                }

                for (NodeIt v(instance.g); v != INVALID; ++v) {
                    for (int scenarioId = 0; scenarioId < instance.nScenarios;
                         scenarioId++) {
                        interiorY[v][scenarioId] = candY[v][scenarioId];
                        UB += (instance.getEdgeRecourseCost(v) /
                               static_cast<double>(instance.nScenarios)) *
                              candY[v][scenarioId];
                    }
                }
            }
            // MIP separation failed, so we stop.
            else {
                break;
            }
        }
        // Otherwise, solve LP again.
        else {
            if (it % 50 == 0) {
                cleanupCuts();
            }

            solveLPAndUpdateSolutions(xValue, yValue, candX, candY);
        }
    }
}

int NodeScenariosModel::addSeparatedCuts(EdgeValueMap &xValue,
                                         NodeVectorValueMap &yValue,
                                         EdgeValueMap &candX,
                                         NodeVectorValueMap &candY,
                                         std::vector<CutData> &separatedCuts) {
    int addedCuts = 0;

    std::vector<size_t> cutOrder(separatedCuts.size());
    std::iota(cutOrder.begin(), cutOrder.end(), 0);
    std::sort(cutOrder.begin(), cutOrder.end(), [&](size_t i, size_t j) {
        return separatedCuts[i].violation > separatedCuts[j].violation;
    });
    int limit = std::min(static_cast<int>(separatedCuts.size()), 5);
    for (int i = 0; i < limit; i++) {
        CutData &cutData = separatedCuts[cutOrder[i]];
        double sign = (cutData.sense == GRB_LESS_EQUAL) ? 1 : -1;
        assert(sign * (cutData.LHS - cutData.RHS) >= 1e-6);

        double LHS_mid = 0.0;
        double LHS_in = 0.0;
        double LHS_out = 0.0;
        GRBLinExpr expr = 0.0;
        for (auto [e, coef] : cutData.edgePairs) {
            assert(e != INVALID);
            assert(std::isfinite(coef));
            expr += coef * x[e];
            LHS_in += coef * interiorX[e];
            LHS_mid += coef * candX[e];
            LHS_out += coef * xValue[e];
        }

        for (auto [scenarioId, scenarioCoef] : cutData.scenarioPairs) {
            for (auto [nodeId, nodeCoef] : cutData.nodePairs) {
                assert(nodeId != instance.depot);
                assert(std::isfinite(nodeCoef));
                Node v = instance.g.nodeFromId(nodeId);
                assert(v != INVALID);
                expr += nodeCoef * y[v][scenarioId];
                LHS_in += nodeCoef * interiorY[v][scenarioId];
                LHS_mid += nodeCoef * candY[v][scenarioId];
                LHS_out += nodeCoef * yValue[v][scenarioId];
            }
        }

        assert(sign * (LHS_in - cutData.RHS) <= 1e-6);
        assert(sign * (LHS_mid - cutData.RHS) >= 1e-6);
        assert(sign * (LHS_out - cutData.RHS) >= 1e-6);

        GRBConstr addedConstraint;
        if (cutPool != nullptr) {
            int cutId = cutPool->addCut(cutData);
            addedConstraint =
                model.addConstr(expr, cutData.sense, cutData.RHS,
                                cutData.name + "?" + std::to_string(cutId));
        } else {
            addedConstraint = model.addConstr(expr, cutData.sense, cutData.RHS);
        }
        addedCuts++;
    }

    return addedCuts;
}

double NodeScenariosModel::buildLagrangianCutAndGetDualObj() {
    // First we use the duals to update node and edge coefs.
    NodeVectorValueMap nodeScenarioCoefs(
        instance.g, std::vector<double>(instance.nScenarios, 0.0));
    EdgeValueMap edgeCoefs(instance.g, 0.0);
    double RHS = 0.0;

    GRBConstr *modelConstraints = model.getConstrs();
    double dualObj = 0.0;
    int nCons = model.get(GRB_IntAttr_NumConstrs);
    for (int i = 0; i < nCons; ++i) {
        std::string name = modelConstraints[i].get(GRB_StringAttr_ConstrName);
        double dual = modelConstraints[i].get(GRB_DoubleAttr_Pi);
        dualObj += dual * modelConstraints[i].get(GRB_DoubleAttr_RHS);

        if (std::abs(dual) <= 1e-6 ||
            name.find("incidence") != std::string::npos ||
            name.find("edgebound") != std::string::npos) {
            continue;
        }

        // Node bound constraint.
        if (name.find("nodebound") != std::string::npos) {
            int pos = name.find_first_of('_');
            std::string tmp = name.substr(pos + 1);
            pos = tmp.find_first_of('_');
            int nodeId = std::stoi(tmp.substr(0, pos));
            int scenarioId = std::stoi(tmp.substr(pos + 1));
            Node v = instance.g.nodeFromId(nodeId);
            nodeScenarioCoefs[v][scenarioId] += dual;
            RHS += dual * modelConstraints[i].get(GRB_DoubleAttr_RHS);
        }

        // SRI constraint.
        else if (name.find("SRI") != std::string::npos) {
            int pos = name.find_first_of('?');
            int cutId = std::stoi(name.substr(pos + 1));
            CutData &cutData = cutPool->getCut(cutId);
            cutData.dual = dual;
            RHS += dual * modelConstraints[i].get(GRB_DoubleAttr_RHS);

            for (auto &[e, coef] : cutData.edgePairs) {
                assert(e != INVALID);
                assert(std::isfinite(coef));
                edgeCoefs[e] += coef * dual;
            }

            for (auto [scenarioId, scenarioCoef] : cutData.scenarioPairs) {
                for (auto [nodeId, nodeCoef] : cutData.nodePairs) {
                    assert(nodeId != instance.depot);
                    assert(std::isfinite(nodeCoef));
                    Node v = instance.g.nodeFromId(nodeId);
                    assert(v != INVALID);
                    nodeScenarioCoefs[v][scenarioId] += nodeCoef * dual;
                }
            }
        }

        // CVRPSEP
        else if (name.find("CVRPSEP") != std::string::npos) {
            int pos = name.find_first_of('?');
            int cutId = std::stoi(name.substr(pos + 1));
            CutData &cutData = cutPool->getCut(cutId);
            cutData.dual = dual;
        }
    }
    delete[] modelConstraints;

    // The node coefs are given by the maximum over each scenario.
    NodeValueMap nodeCoefs(instance.g, 0.0);
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) == instance.depot) {
            continue;
        }

        for (int scenarioId = 0; scenarioId < instance.nScenarios;
             scenarioId++) {
            nodeCoefs[v] =
                std::max(nodeCoefs[v], nodeScenarioCoefs[v][scenarioId]);
        }
    }

    // Now we construct the Lagrangian cut.
    CutData cutData;
    cutData.RHS = RHS;
    cutData.sense = '>';
    cutData.LHS = RHS - 10;
    cutData.name = "LagrangianCut";

    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) == instance.depot) {
            continue;
        }

        if (std::abs(nodeCoefs[v]) >= 1e-6) {
            cutData.nodePairs.push_back(
                {instance.g.id(v), (nodeCoefs[v] * instance.nScenarios) /
                                       instance.getEdgeRecourseCost(v)});
        }
    }

    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        if (std::abs(edgeCoefs[e]) >= 1e-6) {
            cutData.edgePairs.push_back({e, edgeCoefs[e]});
        }
    }

    cutPool->addCut(cutData);

    return dualObj;
}

void NodeScenariosModel::cleanupCuts() {
    GRBConstr *constraints = model.getConstrs();
    int numConstrs = model.get(GRB_IntAttr_NumConstrs);

    std::vector<GRBConstr> toRemove;
    int removedCuts = 0;

    for (int i = 0; i < numConstrs; ++i) {
        std::string name = constraints[i].get(GRB_StringAttr_ConstrName);

        if (name.find("incidence") == std::string::npos) {

            double slack = constraints[i].get(GRB_DoubleAttr_Slack);

            // If slack is large, the cut is not active
            if (std::abs(slack) > 0.1) {
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