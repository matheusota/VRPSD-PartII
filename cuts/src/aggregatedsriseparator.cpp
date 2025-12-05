#include "aggregatedsriseparator.h"

AggregatedSRISeparator::AggregatedSRISeparator(const SVRPInstance &instance,
                                               const Params &params)
    : instance(instance), params(params), optimalRecourseHelper(instance) {
    // Add some scenarios to the critical scenarios set.
    if (params.basicVRPSD) {
        std::vector<int> allScenarios;
        for (int s = 0; s < instance.nScenarios; s++) {
            allScenarios.push_back(s);
        }
        criticalScenarios.insert(allScenarios);
    }

    std::vector<int> topScenarios;
    for (int s = 0; s < 10; s++) {
        topScenarios.push_back(instance.sortedScenarios[s].second);
    }
    criticalScenarios.insert(topScenarios);
}

void AggregatedSRISeparator::buildSupportGraph(const EdgeValueMap &xValue,
                                               Graph &supportGraph,
                                               EdgeValueMap &supportXValue,
                                               NodeIntMap &componentsMap,
                                               int &nComponents) {
    // Add nodes.
    std::vector<Node> nodeById(instance.n, INVALID);
    for (int i = 0; i < instance.n; i++) {
        nodeById[i] = supportGraph.addNode();
    }

    // Add edges to the support graph.
    for (int i = 1; i < instance.n; i++) {
        Node source = nodeById[i];
        Node originalSource = instance.g.nodeFromId(i);
        assert(source != INVALID && originalSource != INVALID);

        for (int j = i + 1; j < instance.n; j++) {
            Node target = nodeById[j];
            Node originalTarget = instance.g.nodeFromId(j);
            assert(target != INVALID && originalTarget != INVALID);

            Edge originalEdge =
                findEdge(instance.g, originalSource, originalTarget);
            assert(originalEdge != INVALID);
            if (xValue[originalEdge] >= EpsForIntegrality) {
                Edge supportEdge = supportGraph.addEdge(source, target);
                supportXValue[supportEdge] = xValue[originalEdge];
            }
        }
    }

    // Delete the depot.
    supportGraph.erase(nodeById[instance.depot]);

    // Check that ids and edges are ok.
    int countNodes = 0;
    for (NodeIt v(supportGraph); v != INVALID; ++v) {
        assert(v != INVALID && supportGraph.id(v) >= 1 &&
               supportGraph.id(v) < instance.n);
        countNodes++;
    }
    assert(countNodes == instance.n - 1);

    for (EdgeIt e(supportGraph); e != INVALID; ++e) {
        Node source = supportGraph.u(e);
        Node originalSource = instance.g.nodeFromId(supportGraph.id(source));
        Node target = supportGraph.v(e);
        Node originalTarget = instance.g.nodeFromId(supportGraph.id(target));
        Edge originalEdge =
            findEdge(instance.g, originalSource, originalTarget);
        assert(originalEdge != INVALID);
        assert(xValue[originalEdge] >= EpsForIntegrality);
    }

    // Find connected components of the support graph.
    nComponents = connectedComponents(supportGraph, componentsMap);
}

int AggregatedSRISeparator::cvrpsepHeuristicSeparation(
    const EdgeValueMap &xValue, const NodeValueMap &recourseValue,
    std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    // Construct support graph.
    Graph supportGraph;
    EdgeValueMap supportXValue(supportGraph, 0.0);
    NodeIntMap componentsMap(supportGraph, 0);
    int nComponents;
    buildSupportGraph(xValue, supportGraph, supportXValue, componentsMap,
                      nComponents);

    // Compute the number of non-zero edges in each connected component.
    std::vector<int> componentNonZero(nComponents, 0);
    for (EdgeIt e(supportGraph); e != INVALID; ++e) {
        Node u = supportGraph.u(e);
        Node v = supportGraph.v(e);

        assert(componentsMap[u] == componentsMap[v] &&
               supportXValue[e] >= EpsForIntegrality);
        assert(componentsMap[u] >= 0 && componentsMap[u] < nComponents);
        componentNonZero[componentsMap[u]]++;
    }

    // Get the customers in each component and the accumulated demands.
    std::vector<std::vector<int>> componentCustomers(nComponents,
                                                     std::vector<int>());
    for (NodeIt v(supportGraph); v != INVALID; ++v) {
        assert(supportGraph.id(v) >= 1 && supportGraph.id(v) < instance.n);
        componentCustomers[componentsMap[v]].push_back(supportGraph.id(v));
    }

    // Add inequalities for each connected component.
    int nSepCuts = 0;
    for (int i = 0; i < nComponents; i++) {
        if (componentCustomers[i].size() <= 1) {
            continue;
        }

        nSepCuts += addCutFromSet(xValue, recourseValue, componentCustomers[i],
                                  separatedCuts);
    }

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return nSepCuts;
}

int AggregatedSRISeparator::addCutFromSet(const EdgeValueMap &xValue,
                                          const NodeValueMap &yValue,
                                          const std::vector<int> &customers,
                                          std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    int addedCuts = 0;

    if (params.dualSetCut) {
        addedCuts +=
            addCutFromSetWithDual(xValue, yValue, customers, separatedCuts);
    } else {
        addedCuts +=
            addCutFromSetSimple(xValue, yValue, customers, separatedCuts);
    }

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return addedCuts;
}

int AggregatedSRISeparator::addCutFromSetSimple(
    const EdgeValueMap &xValue, const NodeValueMap &yValue,
    const std::vector<int> &customers, std::vector<CutData> &separatedCuts) {
    double totalDemand = 0.0;
    for (const int id : customers) {
        Node v = instance.g.nodeFromId(id);
        assert(id >= 1 && id < instance.n && v != INVALID);
        totalDemand += instance.demand[v];

        // There are some instances with zero recourse cost, we ignore those.
        if (instance.getEdgeRecourseCost(v) <= 1e-6) {
            return 0;
        }
    }
    int nVehicles =
        static_cast<int>(std::ceil(totalDemand / instance.capacity));

    // Get x(U).
    double xInside = 0.0;
    for (size_t i = 0; i < customers.size(); i++) {
        for (size_t j = i + 1; j < customers.size(); j++) {
            Edge e = findEdge(instance.g, instance.g.nodeFromId(customers[i]),
                              instance.g.nodeFromId(customers[j]));
            xInside += xValue[e];
        }
    }

    // Check scenarios that we add to the cut.
    int countScenarios = 0;
    double RHS = 0.0;
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        double scenarioSum = 0.0;

        for (const int id : customers) {
            Node v = instance.g.nodeFromId(id);
            scenarioSum += instance.scenariosMatrix[id][scenarioId];
        }

        int currBound = std::ceil(scenarioSum / instance.capacity);
        if (currBound <= 1) {
            continue;
        }

        // When using capacity constraints, we check if k^s(U) \leq \bar{k}(U).
        if (!params.basicVRPSD && currBound <= nVehicles) {
            continue;
        }

        // Check if k^s(U) + x(U) - |U| \leq 0.
        if (currBound + xInside - static_cast<double>(customers.size()) <=
            EpsForIntegrality) {
            continue;
        }

        RHS += currBound - static_cast<int>(customers.size());
        countScenarios++;
    }

    // Construct cut.
    CutData cutData;
    cutData.RHS = RHS;
    cutData.sense = '>';
    cutData.name = "AggregatedSRI-Simple";
    cutData.LHS = 0.0;

    for (size_t i = 0; i < customers.size(); i++) {
        for (size_t j = i + 1; j < customers.size(); j++) {
            Node u = instance.g.nodeFromId(customers[i]);
            Node v = instance.g.nodeFromId(customers[j]);
            Edge e = findEdge(instance.g, u, v);
            assert(e != INVALID);
            cutData.edgePairs.push_back(std::make_pair(e, -countScenarios));
            cutData.LHS += -countScenarios * xValue[e];
        }
    }

    for (const int id : customers) {
        Node v = instance.g.nodeFromId(id);
        assert(v != INVALID);
        cutData.nodePairs.push_back(
            std::make_pair(id, static_cast<double>(instance.nScenarios) /
                                   instance.getEdgeRecourseCost(v)));
        cutData.LHS +=
            (instance.nScenarios * yValue[v]) / instance.getEdgeRecourseCost(v);
        cutData.customers.push_back(id);
    }

    if (cutData.LHS >= cutData.RHS - EpsForIntegrality) {
        return 0;
    } else {
        cutData.violation = cutData.RHS - cutData.LHS;
        nCuts++;
        separatedCuts.push_back(cutData);
        return 1;
    }
}

int AggregatedSRISeparator::addCutFromSetWithDual(
    const EdgeValueMap &xValue, const NodeValueMap &yValue,
    const std::vector<int> &customers, std::vector<CutData> &separatedCuts) {
    double totalDemand = 0.0;
    for (int id : customers) {
        Node v = instance.g.nodeFromId(id);
        assert(id >= 1 && id < instance.n && v != INVALID);
        totalDemand += instance.demand[v];

        // There are some instances with zero recourse cost, we ignore those.
        if (instance.getEdgeRecourseCost(v) <= 1e-6) {
            return 0;
        }
    }
    int nVehicles =
        static_cast<int>(std::ceil(totalDemand / instance.capacity));

    // Sort customers by their recourse values.
    std::vector<std::pair<double, int>> sortedRecourse;
    for (const int id : customers) {
        assert(id >= 1 && id < instance.n);
        Node v = instance.g.nodeFromId(id);
        sortedRecourse.push_back({instance.getEdgeRecourseCost(v), id});
    }
    std::sort(sortedRecourse.begin(), sortedRecourse.end());

    // Construct cut using dual solutions.
    double lowerBound = 0.0;
    double alpha = 0.0;
    NodeValueMap betaCoefs(instance.g, 0.0);
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        // Get sum of scenario demands.
        double scenarioDemand = 0.0;
        for (const int id : customers) {
            scenarioDemand += instance.scenariosMatrix[id][scenarioId];
        }
        int currBound = std::ceil(scenarioDemand / instance.capacity);

        // No failures for this scenario.
        if (currBound <= 1 || currBound <= nVehicles) {
            continue;
        }

        // Get scenario critical item.
        size_t failIndex = static_cast<size_t>(currBound - nVehicles - 1);
        assert(failIndex >= 0 && failIndex < sortedRecourse.size());
        alpha += sortedRecourse[failIndex].first;

        // Compute lower bound and set beta.
        lowerBound += sortedRecourse[failIndex].first;
        for (size_t i = 0; i < sortedRecourse.size(); i++) {
            int id = sortedRecourse[i].second;
            assert(id >= 1 && id < instance.n);
            Node v = instance.g.nodeFromId(id);
            if (i < failIndex) {
                betaCoefs[v] = std::max(betaCoefs[v], sortedRecourse[i].first);
                lowerBound += sortedRecourse[i].first;
            } else {
                betaCoefs[v] =
                    std::max(betaCoefs[v], sortedRecourse[failIndex].first);
            }
        }
    }

    if (lowerBound <= EpsForIntegrality) {
        return 0;
    }

    // Construct cut.
    CutData cutData;
    cutData.RHS =
        lowerBound + alpha * (nVehicles - static_cast<int>(customers.size()));
    cutData.sense = '>';
    cutData.name = "AggregatedSRI-Dual";
    cutData.LHS = 0.0;

    for (size_t i = 0; i < customers.size(); i++) {
        for (size_t j = i + 1; j < customers.size(); j++) {
            Node u = instance.g.nodeFromId(customers[i]);
            Node v = instance.g.nodeFromId(customers[j]);
            Edge e = findEdge(instance.g, u, v);
            assert(e != INVALID);
            cutData.edgePairs.push_back(std::make_pair(e, -alpha));
            cutData.LHS += -alpha * xValue[e];
        }
    }

    for (const int id : customers) {
        Node v = instance.g.nodeFromId(id);
        assert(v != INVALID);
        cutData.nodePairs.push_back(
            std::make_pair(id, (betaCoefs[v] * instance.nScenarios) /
                                   instance.getEdgeRecourseCost(v)));
        cutData.LHS += ((betaCoefs[v] * instance.nScenarios) /
                        instance.getEdgeRecourseCost(v)) *
                       yValue[v];
        cutData.customers.push_back(id);
    }

    if (cutData.LHS >= cutData.RHS - EpsForIntegrality) {
        return 0;
    } else {
        cutData.violation = cutData.RHS - cutData.LHS;
        nCuts++;
        separatedCuts.push_back(cutData);
        return 1;
    }
}

int AggregatedSRISeparator::flowHeuristicSeparation(
    const EdgeValueMap &xValue, const NodeValueMap &recourseValue,
    std::vector<CutData> &separatedCuts) {
    int addedCuts = 0;
    auto started = chrono::high_resolution_clock::now();

    for (auto &scenarios : criticalScenarios) {
        addedCuts += flowSeparationFromScenarios(xValue, recourseValue,
                                                 scenarios, separatedCuts);
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

int AggregatedSRISeparator::flowSeparationFromScenarios(
    const EdgeValueMap &xValue, const NodeValueMap &recourseValue,
    const std::vector<int> &scenarios, std::vector<CutData> &separatedCuts) {
    // Creates a copy of the graph and add a source to it.
    DNode depot = instance.d_g.nodeFromId(instance.depot);
    const double roundFactor = 1e3;
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
        if (v != depot) {
            flowGraph.addArc(newSource, nodeMap[v]);
        }
    }

    // Convert from recourse values to "average aggregated" ys.
    // Also get total demand.
    NodeValueMap y(instance.g, 0.0);
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) != instance.depot) {
            double nodeRecourse = instance.getEdgeRecourseCost(v);
            y[v] = nodeRecourse <= 1e-6
                       ? instance.nScenarios
                       : (recourseValue[v] / nodeRecourse) *
                             static_cast<double>(instance.nScenarios);
            y[v] /= static_cast<double>(scenarios.size());
        }
    }

    // Get top demands.
    double totalDemand = 0.0;
    NodeValueMap topDemands(instance.g, 0.0);
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) != instance.depot) {
            for (int scenarioId : scenarios) {
                topDemands[v] +=
                    instance.scenariosMatrix[instance.g.id(v)][scenarioId];
            }
            topDemands[v] /= static_cast<double>(scenarios.size());
            totalDemand += topDemands[v];
        }
    }

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
            capacityMap[a] =
                static_cast<int>(std::round(topDemands[v2] * roundFactor));
        } else {
            assert(nodeCrossMap[u] != INVALID && nodeCrossMap[v] != INVALID);
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
                capacityValue += instance.capacity * y[u2];
            }
            capacityMap[a] =
                static_cast<int>(std::round(capacityValue * roundFactor));
        }
    }

    // Run max-flow.
    Preflow<Digraph, ArcIntMap> preflow(flowGraph, capacityMap, newSource,
                                        nodeMap[depot]);
    preflow.run();

    // Found violated cut, so we get customers in the cut.
    if (preflow.flowValue() <= totalDemand * roundFactor - 10) {
        std::vector<int> customers;

        double setDemand = 0.0;
        double setY = 0.0;
        for (DNodeIt v(instance.d_g); v != INVALID; ++v) {
            if (preflow.minCut(nodeMap[v])) {
                int id = instance.d_g.id(v);
                assert(id >= 1 && id < instance.n);
                customers.push_back(id);
                setDemand += instance.demand[instance.g.nodeFromId(id)];
                setY += y[instance.g.nodeFromId(id)];
            }
        }

        int addedCuts =
            addCutFromSet(xValue, recourseValue, customers, separatedCuts);
        // assert(addedCuts >= 1);
        return addedCuts;
    }

    return 0;
}

int AggregatedSRISeparator::mipSeparation(const EdgeValueMap &xValue,
                                          const NodeValueMap &recourseValue,
                                          std::vector<CutData> &separatedCuts) {
    GRBEnv env(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    GRBModel model(env);

    model.set(GRB_StringAttr_ModelName, "GurobiModel");
    model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    model.set(GRB_DoubleParam_TimeLimit, 30);
    model.set(GRB_IntParam_MIPFocus, 1);
    model.set(GRB_IntParam_Threads, 1);
    env.set(GRB_IntParam_LogToConsole, 0);

    // Set Gurobi variables.
    std::vector<GRBVar> scenarioVars(instance.nScenarios);
    std::vector<GRBVar> scenarioVehicles(instance.nScenarios);
    for (int s = 0; s < instance.nScenarios; s++) {
        scenarioVars[s] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                       "scenarioVar_" + std::to_string(s));
        scenarioVehicles[s] =
            model.addVar(-1.0, GRB_INFINITY, 0.0, GRB_INTEGER,
                         "scenarioVehicles_" + std::to_string(s));
    }

    NodeGRBVarMap customersReal(instance.g);
    NodeGRBVarMap customersBinary(instance.g);
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) == instance.depot) {
            continue;
        }

        customersBinary[v] =
            model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                         "customersBinary_" + std::to_string(instance.g.id(v)));
        customersReal[v] =
            model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
                         "customersReal_" + std::to_string(instance.g.id(v)));
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
            model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
                         "edgesReal_" + std::to_string(instance.g.id(u)) + "," +
                             std::to_string(instance.g.id(v)));
    }
    model.update();

    // Objective bound constraint.
    GRBLinExpr objExpr = 0.0;
    for (int s = 0; s < instance.nScenarios; s++) {
        objExpr += 1.0 + scenarioVehicles[s];
    }

    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) == instance.depot) {
            continue;
        }

        // Need to treat zero recourse value.
        double coef = 0.0;
        double recourseCost = instance.getEdgeRecourseCost(v);

        if (recourseCost >= 1e-6) {
            coef = (recourseValue[v] * instance.nScenarios) / recourseCost;
        } else if (recourseCost <= 1e-6 && recourseValue[v] >= 1e-6) {
            model.addConstr(customersBinary[v] == 0);
        }

        objExpr -= coef * customersBinary[v] + customersReal[v];
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

    // Scenario constraints.
    GRBLinExpr scenarioSumExpr = 0.0;
    int nVehicles = params.basicVRPSD ? instance.n - 1 : instance.k;
    for (int s = 0; s < instance.nScenarios; s++) {
        scenarioSumExpr += scenarioVars[s];

        GRBLinExpr loadExpr = -0.1;
        GRBLinExpr avgLoadExpr = 0.0;
        for (NodeIt v(instance.g); v != INVALID; ++v) {
            if (instance.g.id(v) == instance.depot) {
                continue;
            }
            double demand = instance.scenariosMatrix[instance.g.id(v)][s];
            loadExpr += demand * customersBinary[v];
            avgLoadExpr +=
                (instance.demand[v] / instance.capacity) * customersBinary[v];
        }

        model.addConstr(instance.capacity * scenarioVehicles[s] <= loadExpr);
        model.addConstr(scenarioVehicles[s] <=
                        nVehicles * scenarioVars[s] - 1.0);

        // If not using basic, we need to account for the average number of
        // vehicles.
        if (!params.basicVRPSD) {
            model.addConstr(scenarioVehicles[s] >=
                            avgLoadExpr - nVehicles * (1.0 - scenarioVars[s]));
        }
    }

    // Add customer constraints.
    GRBLinExpr sizeExpr = 0.0;
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) == instance.depot) {
            continue;
        }

        sizeExpr += customersBinary[v];
        model.addConstr(customersReal[v] >=
                        scenarioSumExpr -
                            instance.nScenarios * (1.0 - customersBinary[v]));
        model.addConstr(customersReal[v] <=
                        instance.nScenarios * customersBinary[v]);
    }
    model.addConstr(sizeExpr >= 2);

    // Add edge constraints.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        Node u = instance.g.u(e);
        Node v = instance.g.v(e);
        if (instance.g.id(u) == instance.depot ||
            instance.g.id(v) == instance.depot) {
            continue;
        }

        model.addConstr(edgesReal[e] <= customersReal[u]);
        model.addConstr(edgesReal[e] <= customersReal[v]);
        model.addConstr(edgesReal[e] <= scenarioSumExpr);
    }

    model.update();
    model.optimize();

    int addedCuts = 0;
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        std::vector<int> customers;
        for (NodeIt v(instance.g); v != INVALID; ++v) {
            if (instance.g.id(v) == instance.depot) {
                continue;
            }

            if (customersBinary[v].get(GRB_DoubleAttr_X) >= 0.1) {
                customers.push_back(instance.g.id(v));
            }
        }

        int addedCuts = addCutFromSetSimple(xValue, recourseValue, customers,
                                            separatedCuts);
        assert(addedCuts >= 1);
    }

    return addedCuts;
}