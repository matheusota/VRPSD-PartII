#include "sriseparator.h"

SRISeparator::SRISeparator(const SVRPInstance &instance, const Params &params,
                           CVRPSEPSeparator &cvrpsepSeparator)
    : instance(instance), params(params), cvrpsepSeparator(cvrpsepSeparator) {}

int SRISeparator::flowSeparation(const EdgeValueMap &xValue,
                                 const NodeVectorValueMap &yValue,
                                 std::vector<CutData> &separatedCuts) {
    // Creates a copy of the graph and add a source to it.
    DNode depot = instance.d_g.nodeFromId(instance.depot);
    double roundFactor = 1e6;
    Digraph flowGraph;
    DNodeDNodeMap nodeMap(instance.d_g);
    DNodeDNodeMap nodeCrossMap(flowGraph);
    ArcArcMap arcMap(instance.d_g);
    ArcArcMap arcCrossMap(flowGraph);
    digraphCopy(instance.d_g, flowGraph)
        .nodeRef(nodeMap)
        .nodeCrossRef(nodeCrossMap)
        .arcRef(arcMap)
        .arcCrossRef(arcCrossMap)
        .run();

    DNode newSource = flowGraph.addNode();
    for (DNodeIt v(instance.d_g); v != INVALID; ++v) {
        if (instance.d_g.id(v) != instance.depot) {
            flowGraph.addArc(newSource, nodeMap[v]);
        }
    }

    // Run MaxFlow-MinCut for each scenario.
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        // Create arc capacities map.
        ArcIntMap capacityMap(flowGraph, 0);
        for (ArcIt a(flowGraph); a != INVALID; ++a) {
            DNode u = flowGraph.source(a);
            DNode v = flowGraph.target(a);

            if (u == newSource) {
                assert(nodeCrossMap[v] != INVALID);
                int id = instance.d_g.id(nodeCrossMap[v]);
                assert(id >= 1 && id < instance.n);
                Node v2 = instance.g.nodeFromId(id);
                capacityMap[a] = static_cast<int>(std::round(
                    instance.scenariosMatrix[id][scenarioId] * roundFactor));
            } else {
                assert(nodeCrossMap[u] != INVALID &&
                       nodeCrossMap[v] != INVALID);
                int id_u = instance.d_g.id(nodeCrossMap[u]);
                int id_v = instance.d_g.id(nodeCrossMap[v]);
                assert(id_u >= 0 && id_u < instance.n && id_v >= 0 &&
                       id_v < instance.n);
                Node u2 = instance.g.nodeFromId(id_u);
                Node v2 = instance.g.nodeFromId(id_v);
                Edge e = findEdge(instance.g, u2, v2);
                assert(e != INVALID);

                double capacityValue = (0.5 * instance.capacity) * xValue[e];
                if (id_v == instance.depot) {
                    capacityValue += instance.capacity * yValue[u2][scenarioId];
                }
                capacityMap[a] =
                    static_cast<int>(std::round(capacityValue * roundFactor));
            }
        }

        // Run max-flow.
        Preflow<Digraph, ArcIntMap> preflow(flowGraph, capacityMap, newSource,
                                            nodeMap[depot]);
        preflow.run();

        // Get total scenarioDemand.
        double scenarioDemand = 0.0;
        for (int id = 1; id < instance.n; id++) {
            scenarioDemand += instance.scenariosMatrix[id][scenarioId];
        }

        // Found violated cut, so we get customers in the cut.
        if (preflow.flowValue() <= scenarioDemand * roundFactor - 10) {
            std::vector<int> customers;

            for (DNodeIt v(instance.d_g); v != INVALID; ++v) {
                if (preflow.minCut(nodeMap[v])) {
                    int id = instance.d_g.id(v);
                    assert(id >= 1 && id < instance.n);
                    customers.push_back(id);
                }
            }

            int addedCuts =
                addSRIFromCustomerSet(xValue, yValue, customers, separatedCuts);

            assert(addedCuts >= 1);
            return addedCuts;
        }
    }

    return 0;
}

int SRISeparator::heuristicSeparation(const EdgeValueMap &xValue,
                                      const NodeVectorValueMap &yValue,
                                      std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    int addedCuts = 0;
    for (int i = 0; i < instance.nScenarios; i++) {
        int scenarioId = instance.sortedScenarios[i].second;

        // Call CVRPSEP for each scenario to update customer sets.
        if (cvrpsepSeparator.updateCustomerSetsForScenario(xValue,
                                                           scenarioId)) {
            const std::vector<std::vector<int>> &customerSets =
                cvrpsepSeparator.getCustomerSets();

            // Try to separate node recourse cuts for the customer sets.
            for (size_t i = 0; i < customerSets.size(); i++) {
                addedCuts += addSRIFromCustomerSet(
                    xValue, yValue, customerSets[i], separatedCuts);
            }

            if (addedCuts > 0) {
                break;
            }
        }
    }

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return addedCuts;
}

int SRISeparator::mipSeparation(const EdgeValueMap &xValue,
                                const NodeVectorValueMap &yValue,
                                std::vector<CutData> &separatedCuts,
                                SRI_SEPARATION_STATUS &status) {
    auto started = chrono::high_resolution_clock::now();

    int addedCuts = 0;
    for (int i = 0; i < instance.nScenarios; i++) {
        int scenarioId = instance.sortedScenarios[i].second;

        std::vector<int> customers;
        SRI_SEPARATION_STATUS sepStatus =
            singleScenarioMIPModel(xValue, yValue, scenarioId, customers);
        switch (sepStatus) {
        case SRI_SEPARATION_STATUS::NO_CUTS: {
            break;
        }
        case SRI_SEPARATION_STATUS::FOUND_CUTS: {
            // Try to separate node recourse cuts for the customer sets.
            addedCuts +=
                addSRIFromCustomerSet(xValue, yValue, customers, separatedCuts);
            assert(addedCuts >= 1);
            status = SRI_SEPARATION_STATUS::FOUND_CUTS;
            break;
        }
        case SRI_SEPARATION_STATUS::FAILED: {
            if (status == SRI_SEPARATION_STATUS::NO_CUTS) {
                status = SRI_SEPARATION_STATUS::FAILED;
            }
        }
        }

        if (addedCuts >= 1) {
            break;
        }
    }

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return addedCuts;
}

int SRISeparator::separateSRIFromCustomerSet(
    const EdgeValueMap &xValue, const NodeVectorValueMap &yValue,
    const std::vector<int> &customers, std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    int addedCuts =
        addSRIFromCustomerSet(xValue, yValue, customers, separatedCuts);

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return addedCuts;
}

int SRISeparator::integerSeparation(const EdgeValueMap &xValue,
                                    const NodeVectorValueMap &yValue,
                                    std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    // This should only be called for integer solutions.
    EdgeIntMap xIntMap(instance.g);

    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        assert(std::abs(xValue[e] - std::round(xValue[e])) <=
               EpsForIntegrality);
        xIntMap[e] = static_cast<int>(std::round(xValue[e]));
    }

    // Get routes.
    SVRPSolution solution;
    solution.setRoutesFromSolution(instance, xIntMap);

    // Call separation for each route.
    int addedCuts = 0;
    for (const auto &route : solution.routes) {
        addedCuts += addSRICutsFromRoute(xValue, yValue, route, separatedCuts);
    }

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return addedCuts;
}

int SRISeparator::addSRIFromCustomerSet(const EdgeValueMap &xValue,
                                        const NodeVectorValueMap &yValue,
                                        const std::vector<int> &customers,
                                        std::vector<CutData> &separatedCuts) {
    // Get x(U).
    double xInside = 0.0;
    for (size_t i = 0; i < customers.size(); i++) {
        for (size_t j = i + 1; j < customers.size(); j++) {
            Edge e = findEdge(instance.g, instance.g.nodeFromId(customers[i]),
                              instance.g.nodeFromId(customers[j]));
            xInside += xValue[e];
        }
    }

    // Get scenario bounds.
    std::vector<std::pair<int, int>> scenarioBounds;
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        double scenarioSum = 0.0;
        double customersY = 0.0;

        for (const int id : customers) {
            Node v = instance.g.nodeFromId(id);
            scenarioSum += instance.scenariosMatrix[id][scenarioId];
            customersY += yValue[v][scenarioId];
        }

        int currBound = std::ceil(scenarioSum / instance.capacity);
        // Check if y^s(U) < k^s(U) + x(U) - |U|
        if (customersY <= currBound + xInside -
                              static_cast<double>(customers.size()) -
                              EpsForIntegrality) {
            scenarioBounds.push_back({scenarioId, currBound});
        }
    }

    if (scenarioBounds.empty()) {
        return 0;
    }

    return addCut(xValue, yValue, scenarioBounds, customers, separatedCuts);
}

// Add the inequality y(U) \geq \sum_{s \in S'} (k_s(U) - x(\delta(U)) / 2).
// In fact, since 2 * x(E(U)) + x(\delta(U)) = 2 * |U|, we add the
// inequality y(U) \geq \sum_{s \in S'} (k_s(U) - |U| + x(E(U))) instead.
int SRISeparator::addCut(const EdgeValueMap &xValue,
                         const NodeVectorValueMap &yValue,
                         const std::vector<std::pair<int, int>> &scenarioBounds,
                         const std::vector<int> &customers,
                         std::vector<CutData> &separatedCuts) {
    separatedCuts.emplace_back(CutData());
    CutData &cutData = separatedCuts.back();
    cutData.RHS = 0.0;
    cutData.sense = '>';
    cutData.name = "SRI";
    cutData.LHS = 0.0;
    int nScenariosCut = static_cast<int>(scenarioBounds.size());

    for (auto &[scenarioId, bound] : scenarioBounds) {
        assert(scenarioId >= 0 && scenarioId < instance.nScenarios);
        cutData.RHS += bound - static_cast<int>(customers.size());
        cutData.scenarioPairs.push_back(std::make_pair(scenarioId, 1.0));

        for (const int id : customers) {
            assert(id >= 1 && id < instance.n);
            Node v = instance.g.nodeFromId(id);
            cutData.LHS += yValue[v][scenarioId];
        }
    }

    for (size_t i = 0; i < customers.size(); i++) {
        for (size_t j = i + 1; j < customers.size(); j++) {
            Node u = instance.g.nodeFromId(customers[i]);
            Node v = instance.g.nodeFromId(customers[j]);
            Edge e = findEdge(instance.g, u, v);
            assert(e != INVALID);
            cutData.edgePairs.push_back(std::make_pair(e, -nScenariosCut));
            cutData.LHS += -nScenariosCut * xValue[e];
        }
    }

    for (const int id : customers) {
        assert(id >= 1 && id < instance.n);
        cutData.nodePairs.push_back(std::make_pair(id, 1.0));
        cutData.customers.push_back(id);
    }

    if (cutData.LHS >= cutData.RHS - EpsForIntegrality) {
        separatedCuts.pop_back();
        return 0;
    } else {
        cutData.violation = cutData.RHS - cutData.LHS;
        nCuts++;
        return 1;
    }
}

SRI_SEPARATION_STATUS SRISeparator::singleScenarioMIPModel(
    const EdgeValueMap &xValue, const NodeVectorValueMap &yValue,
    int scenarioId, std::vector<int> &customers) {
    GRBEnv env(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    GRBModel model(env);

    model.set(GRB_StringAttr_ModelName, "GurobiModel");
    model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    model.set(GRB_DoubleParam_TimeLimit, 30);
    model.set(GRB_IntParam_MIPFocus, 1);
    // model.set(GRB_DoubleParam_MIPGap, 1e-3);
    model.set(GRB_IntParam_Threads, 1);
    env.set(GRB_IntParam_LogToConsole, 0);

    // Set Gurobi variables.
    GRBVar scenarioVehicles =
        model.addVar(-1.0, GRB_INFINITY, 0.0, GRB_INTEGER,
                     "scenarioVehicles_" + std::to_string(scenarioId));

    NodeGRBVarMap customersBinary(instance.g);
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) == instance.depot) {
            continue;
        }

        customersBinary[v] =
            model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                         "customersBinary_" + std::to_string(instance.g.id(v)));
    }

    EdgeGRBVarMap edgesReal(instance.g);
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        Node u = instance.g.u(e);
        Node v = instance.g.v(e);
        if (instance.g.id(u) == instance.depot ||
            instance.g.id(v) == instance.depot) {
            continue;
        }

        edgesReal[e] =
            model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS,
                         "edgesReal_" + std::to_string(instance.g.id(u)) + "," +
                             std::to_string(instance.g.id(v)));
    }
    model.update();

    // Bound constraint.
    GRBLinExpr objExpr = 1 + scenarioVehicles;

    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) == instance.depot) {
            continue;
        }

        double coef = yValue[v][scenarioId];

        objExpr -= (coef + 1) * customersBinary[v];
    }

    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        Node u = instance.g.u(e);
        Node v = instance.g.v(e);
        if (instance.g.id(u) == instance.depot ||
            instance.g.id(v) == instance.depot) {
            continue;
        }

        objExpr += xValue[e] * edgesReal[e];
    }

    model.addConstr(objExpr >= 0.1);

    // Add scenario constraints.
    GRBLinExpr loadExpr = -0.1;
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) == instance.depot) {
            continue;
        }
        double demand = instance.scenariosMatrix[instance.g.id(v)][scenarioId];
        loadExpr += demand * customersBinary[v];
    }

    model.addConstr(instance.capacity * scenarioVehicles <= loadExpr);

    // Add customer constraints.
    GRBLinExpr customersExpr = 0;
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) == instance.depot) {
            continue;
        }
        customersExpr += customersBinary[v];
    }
    model.addConstr(customersExpr >= 2);

    // Add edge constraints.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        Node u = instance.g.u(e);
        Node v = instance.g.v(e);
        if (instance.g.id(u) == instance.depot ||
            instance.g.id(v) == instance.depot) {
            continue;
        }

        model.addConstr(edgesReal[e] <= customersBinary[u]);
        model.addConstr(edgesReal[e] <= customersBinary[v]);
        model.addConstr(edgesReal[e] >=
                        customersBinary[u] + customersBinary[v] - 1.0);
    }

    model.update();
    model.optimize();

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        for (NodeIt v(instance.g); v != INVALID; ++v) {
            if (instance.g.id(v) == instance.depot) {
                continue;
            }

            if (customersBinary[v].get(GRB_DoubleAttr_X) >= 0.1) {
                customers.push_back(instance.g.id(v));
            }
        }

        return SRI_SEPARATION_STATUS::FOUND_CUTS;
    } else if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
        return SRI_SEPARATION_STATUS::NO_CUTS;
    } else {
        return SRI_SEPARATION_STATUS::FAILED;
    }
}

int SRISeparator::addSRICutsFromRoute(const EdgeValueMap &xValue,
                                      const NodeVectorValueMap &yValue,
                                      const std::vector<int> &route,
                                      std::vector<CutData> &separatedCuts) {
    int addedCuts = 0;

    // Loop through all subroutes.
    for (size_t i = 0; i < route.size(); i++) {
        std::vector<int> customers;
        std::vector<int> accDemand(instance.nScenarios, 0);
        std::vector<double> ySum(instance.nScenarios, 0.0);

        for (size_t j = i; j < route.size(); j++) {
            customers.push_back(route[j]);

            // Get scenario bounds for this subroute.
            std::vector<std::pair<int, int>> scenarioBounds;
            for (int scenarioId = 0; scenarioId < instance.nScenarios;
                 scenarioId++) {
                assert(route[j] >= 0 && route[j] < instance.n);
                accDemand[scenarioId] +=
                    instance.scenariosMatrix[route[j]][scenarioId];
                ySum[scenarioId] +=
                    yValue[instance.g.nodeFromId(route[j])][scenarioId];

                int currBound =
                    std::ceil(static_cast<double>(accDemand[scenarioId]) /
                              instance.capacity);
                double LHS = ySum[scenarioId] + 1.0;
                if (currBound >= LHS + EpsForIntegrality) {
                    scenarioBounds.push_back(
                        std::make_pair(scenarioId, currBound));
                }
            }

            if (scenarioBounds.empty()) {
                continue;
            }

            addedCuts += addCut(xValue, yValue, scenarioBounds, customers,
                                separatedCuts);

            if (addedCuts >= 5) {
                break;
            }
        }
    }

    return addedCuts;
}

int SRISeparator::integerSeparation2(const EdgeValueMap &xValue,
                                     const NodeValueMap &yValue,
                                     std::vector<CutData> &separatedCuts) {
    // This should only be called for integer solutions.
    EdgeIntMap xIntMap(instance.g);

    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        assert(std::abs(xValue[e] - std::round(xValue[e])) <=
               EpsForIntegrality);
        xIntMap[e] = static_cast<int>(std::round(xValue[e]));
    }

    // Get routes.
    SVRPSolution solution;
    solution.setRoutesFromSolution(instance, xIntMap);

    // Call separation for each route.
    int addedCuts = 0;
    for (const auto &route : solution.routes) {
        double recourseCost =
            instance.optimalRecourseHelper.getRouteRecourseCost(route);
        double LHS = 0.0;
        for (int id : route) {
            Node v = instance.g.nodeFromId(id);
            LHS += instance.getEdgeRecourseCost(v) * yValue[v];
        }

        if (LHS <= recourseCost - EpsForIntegrality) {
            addedCuts += addBendersCut(xValue, yValue, route, separatedCuts);
            if (addedCuts > 0) {
                break;
            }
        }
    }

    return addedCuts;
}

int SRISeparator::addBendersCut(const EdgeValueMap &xValue,
                                const NodeValueMap &yValue,
                                const std::vector<int> &route,
                                std::vector<CutData> &separatedCuts) {
    // Get tight subroutes for each scenario.
    std::vector<std::vector<std::unordered_set<int>>> tightSubroutes;
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        std::vector<std::unordered_set<int>> subroutes;
        instance.optimalRecourseHelper.getRouteRecourseCostInScenarioWithLP(
            route, scenarioId, subroutes);
        tightSubroutes.push_back(subroutes);
    }

    GRBEnv env(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    GRBModel model(env);

    model.set(GRB_StringAttr_ModelName, "GurobiModel");
    model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    model.set(GRB_DoubleParam_TimeLimit, 300);
    model.set(GRB_IntParam_Threads, 1);
    env.set(GRB_IntParam_LogToConsole, 0);

    // Create alpha variables for each scenario.
    std::vector<std::vector<GRBVar>> alpha;

    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        alpha.push_back(std::vector<GRBVar>(tightSubroutes[scenarioId].size()));
        for (size_t r = 0; r < tightSubroutes[scenarioId].size(); r++) {
            alpha[scenarioId][r] =
                model.addVar(0.0, 1.0, 1.0, GRB_CONTINUOUS,
                             "alpha_" + std::to_string(scenarioId) + "_" +
                                 std::to_string(r));
        }
    }

    // Create beta variables for each customer.
    std::vector<GRBVar> beta(route.size());

    for (size_t i = 0; i < route.size(); i++) {
        Node v = instance.g.nodeFromId(route[i]);
        beta[i] = model.addVar(0.0, GRB_INFINITY, -yValue[v], GRB_CONTINUOUS,
                               "beta_" + std::to_string(route[i]));
    }
    model.update();

    // Add subroute constraints.
    for (size_t i = 0; i < route.size(); i++) {
        for (int scenarioId = 0; scenarioId < instance.nScenarios;
             scenarioId++) {
            GRBLinExpr expr = 0.0;

            for (size_t r = 0; r < tightSubroutes[scenarioId].size(); r++) {
                if (tightSubroutes[scenarioId][r].find(route[i]) !=
                    tightSubroutes[scenarioId][r].end()) {
                    expr += alpha[scenarioId][r];
                }
            }

            model.addConstr(expr <= beta[i]);
        }
    }

    model.update();
    model.optimize();

    if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
        return 0;
    }

    double objVal = model.get(GRB_DoubleAttr_ObjVal);
    if (objVal <= EpsForIntegrality) {
        std::cout << "No Benders cut to separate." << std::endl;
        return 0;
    }

    // Construct cut.
    EdgeValueMap edgeCoefs(instance.g, 0.0);
    NodeValueMap nodeCoefs(instance.g, 0.0);
    double RHS = 0.0;
    double tmp = 0.0;

    for (size_t i = 0; i < route.size(); i++) {
        Node v = instance.g.nodeFromId(route[i]);
        nodeCoefs[v] = beta[i].get(GRB_DoubleAttr_X);
    }

    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        for (size_t r = 0; r < tightSubroutes[scenarioId].size(); r++) {
            double currAlpha = alpha[scenarioId][r].get(GRB_DoubleAttr_X);

            if (currAlpha >= EpsForIntegrality) {
                RHS += currAlpha *
                       (2.0 - static_cast<double>(
                                  tightSubroutes[scenarioId][r].size()));
                tmp -=
                    static_cast<double>(tightSubroutes[scenarioId][r].size());

                for (int id1 : tightSubroutes[scenarioId][r]) {
                    for (int id2 : tightSubroutes[scenarioId][r]) {
                        if (id1 != id2) {
                            Node u = instance.g.nodeFromId(id1);
                            Node v = instance.g.nodeFromId(id2);
                            Edge e = findEdge(instance.g, u, v);
                            edgeCoefs[e] -= 0.5 * currAlpha;
                            tmp += 0.5 * currAlpha * xValue[e];
                        }
                    }
                }
            }
        }
    }

    separatedCuts.emplace_back(CutData());
    CutData &cutData = separatedCuts.back();
    cutData.RHS = RHS;
    cutData.sense = '>';
    cutData.name = "BendersCut";
    cutData.LHS = 0.0;

    for (const int id : route) {
        assert(id >= 1 && id < instance.n);
        Node v = instance.g.nodeFromId(id);
        cutData.LHS += nodeCoefs[v] * yValue[v];
        cutData.nodePairs.push_back(std::make_pair(id, nodeCoefs[v]));
        cutData.customers.push_back(id);
    }

    for (size_t i = 0; i < route.size() - 1; i++) {
        Node u = instance.g.nodeFromId(route[i]);
        Node v = instance.g.nodeFromId(route[i + 1]);
        Edge e = findEdge(instance.g, u, v);
        assert(e != INVALID);
        cutData.edgePairs.push_back(std::make_pair(e, edgeCoefs[e]));
        cutData.LHS += edgeCoefs[e] * xValue[e];
    }

    if (cutData.LHS >= cutData.RHS - EpsForIntegrality) {
        separatedCuts.pop_back();
        return 0;
    } else {
        cutData.violation = cutData.RHS - cutData.LHS;
        nCuts++;
        return 1;
    }
}