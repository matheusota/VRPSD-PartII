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

    // Create the digraph for the network flow heuristic.
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        flowGraph.addNode();
    }
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        Node u = instance.g.u(e);
        Node v = instance.g.v(e);
        DNode du = flowGraph.nodeFromId(instance.g.id(u));
        DNode dv = flowGraph.nodeFromId(instance.g.id(v));
        flowGraph.addArc(du, dv);
        flowGraph.addArc(dv, du);
    }

    // Artificial source and arcs.
    newSource = flowGraph.addNode();
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) != instance.depot) {
            DNode dv = flowGraph.nodeFromId(instance.g.id(v));
            flowGraph.addArc(newSource, dv);
        }
    }
}

int AggregatedSRISeparator::addCutFromSet(const EdgeValueMap &xValue,
                                          const NodeValueMap &yValue,
                                          const std::vector<int> &customers,
                                          std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    int addedCuts =
        addCutFromSetSimple(xValue, yValue, customers, separatedCuts);
    // int addedCuts =
    //     addCutFromSetWithDual(xValue, yValue, customers, separatedCuts);

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

        RHS += currBound - static_cast<double>(customers.size());
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
            Edge e = findEdge(instance.g, instance.g.nodeFromId(customers[i]),
                              instance.g.nodeFromId(customers[j]));
            assert(e != INVALID);
            cutData.edgePairs.push_back({e, -countScenarios});
            cutData.LHS += -countScenarios * xValue[e];
        }
    }

    for (const int id : customers) {
        Node v = instance.g.nodeFromId(id);
        assert(v != INVALID);
        cutData.nodePairs.push_back({id, 1.0});
        cutData.LHS += yValue[v];
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
    const EdgeValueMap &xValue, const NodeValueMap &yValue,
    std::vector<CutData> &separatedCuts) {
    int addedCuts = 0;
    auto started = chrono::high_resolution_clock::now();

    for (auto &scenarios : criticalScenarios) {
        addedCuts += flowSeparationFromScenarios(xValue, yValue, scenarios,
                                                 separatedCuts);
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
    const EdgeValueMap &xValue, const NodeValueMap &yValue,
    const std::vector<int> &scenarios, std::vector<CutData> &separatedCuts) {
    NodeValueMap scaledY(instance.g, 0.0);
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) != instance.depot) {
            if (instance.getEdgeRecourseCost(v) <= 1e-6) {
                scaledY[v] = (instance.n - 1) * instance.nScenarios;
            } else {
                // We divide here because we take the average of the scenario
                // demands.
                scaledY[v] = yValue[v] / static_cast<double>(scenarios.size());
            }
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
        DNode du = flowGraph.source(a);
        DNode dv = flowGraph.target(a);

        if (du == newSource) {
            int id = flowGraph.id(dv);
            assert(id >= 1 && id < instance.n);
            Node v = instance.g.nodeFromId(id);
            capacityMap[a] =
                static_cast<int>(std::ceil(topDemands[v] * flowRoundFactor));
        } else {
            int idu = flowGraph.id(du);
            int idv = flowGraph.id(dv);
            assert(idu >= 0 && idu < instance.n && idv >= 0 &&
                   idv < instance.n);
            Node u = instance.g.nodeFromId(idu);
            Node v = instance.g.nodeFromId(idv);
            Edge e = findEdge(instance.g, u, v);
            assert(e != INVALID);

            double capacityValue = (0.5 * instance.capacity) * xValue[e];
            if (idv == instance.depot) {
                capacityValue += instance.capacity * scaledY[u];
            }
            capacityMap[a] =
                static_cast<int>(std::ceil(capacityValue * flowRoundFactor));
        }
    }

    // Run max-flow.
    Preflow<Digraph, ArcIntMap> preflow(flowGraph, capacityMap, newSource,
                                        flowGraph.nodeFromId(instance.depot));
    preflow.run();

    // Found violated cut, so we get customers in the cut.
    int addedCuts = 0;
    if (preflow.flowValue() <= totalDemand * flowRoundFactor - 10) {
        std::vector<int> customers;
        double customersY = 0.0;
        double customersDemand = 0.0;
        double xInside = 0.0;

        for (DNodeIt v(flowGraph); v != INVALID; ++v) {
            if (preflow.minCut(v) && v != newSource) {
                int id = flowGraph.id(v);
                assert(id >= 1 && id < instance.n);
                customers.push_back(id);
                Node v2 = instance.g.nodeFromId(id);
                customersY += yValue[v2];
                customersDemand += topDemands[v2];
            }
        }

        assert(!customers.empty());
        for (size_t i = 0; i < customers.size(); i++) {
            for (size_t j = i + 1; j < customers.size(); j++) {
                Node u = instance.g.nodeFromId(customers[i]);
                Node v = instance.g.nodeFromId(customers[j]);
                Edge e = findEdge(instance.g, u, v);
                xInside += xValue[e];
            }
        }

        assert(customersY <=
               static_cast<double>(scenarios.size()) *
                       ((customersDemand / instance.capacity) + xInside -
                        static_cast<double>(customers.size())) -
                   EpsForIntegrality);

        addedCuts =
            addCutFromSetSimple(xValue, yValue, customers, separatedCuts);
    }

    return addedCuts;
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
        if (currBound <= 1) {
            continue;
        }

        // When using capacity constraints, we check if k^s(U) \leq \bar{k}(U).
        if (!params.basicVRPSD && currBound <= nVehicles) {
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