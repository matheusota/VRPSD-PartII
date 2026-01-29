#include "partialrouteseparator.h"
#include "partialroutesbuilder.h"

PartialRouteSeparator::PartialRouteSeparator(SVRPInstance &instance,
                                             const Params &params)
    : instance(instance), params(params) {}

int PartialRouteSeparator::separatePartialRouteCuts(
    const EdgeValueMap &xValue, const NodeValueMap &recourseValue,
    const std::vector<PartialRoute> &partialRoutes,
    std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();
    int nSepCuts = 0;

    // Add inequalities for each partial route.
    for (const PartialRoute &partialRoute : partialRoutes) {
        assert(partialRoute.excess >= EpsForIntegrality);

        // Get partial route costs.
        NodeValueMap nodeRecourseMap(instance.g, 0.0);
        double partialRouteClassicalCost =
            instance.classicalRecourseHelper
                .getPartialRouteRecourseAndFailureNodes(partialRoute,
                                                        nodeRecourseMap);
        std::unordered_map<std::pair<int, int>, double,
                           boost::hash<std::pair<int, int>>>
            betaDuals;
        std::unordered_map<std::tuple<int, int, int>, double,
                           boost::hash<std::tuple<int, int, int>>>
            alphaDuals;
        double partialRouteOptimalCost = 0.0;
        if (params.scenarioOptimalPRCuts) {
            partialRouteOptimalCost =
                instance.optimalRecourseHelper.getPartialRouteRecourseCost(
                    partialRoute, betaDuals, alphaDuals);
        }

        if (partialRouteOptimalCost <= EpsForIntegrality &&
            partialRouteClassicalCost <= EpsForIntegrality) {
            continue;
        }

        // Pass through all nodes in the partial route.
        double currCost = 0.0;
        NodeValueMap nodeCoefs(instance.g, 0.0);
        for (const PartialRouteEntry &entry : partialRoute.entries) {
            for (int id : entry.vertices) {
                Node v = instance.g.nodeFromId(id);
                currCost += recourseValue[v];
                nodeCoefs[v] = 1.0;
            }
        }

        // Compare costs.
        bool useSimpleAdherenceCuts = false;
        double recourseCostToUse = partialRouteClassicalCost;
        if (params.policy == SCENARIO_OPTIMAL ||
            currCost + EpsForIntegrality <= partialRouteOptimalCost) {
            recourseCostToUse = partialRouteOptimalCost;
            useSimpleAdherenceCuts = true;
        }

        if (currCost >= recourseCostToUse - EpsForIntegrality) {
            continue;
        }

        // Setup the coefficients for the inequality.
        double RHS = 0.0;
        EdgeValueMap edgeCoefs(instance.g, 0.0);
        if (params.projectedSRI) {
            setCutCoefficientsFromDual(partialRoute, xValue, recourseCostToUse,
                                       betaDuals, alphaDuals, edgeCoefs,
                                       nodeCoefs, RHS);
        } else {
            setCutCoefficientsNew(partialRoute, recourseCostToUse, edgeCoefs,
                                  RHS, useSimpleAdherenceCuts,
                                  params.scenarioOptimalPRCuts);
        }

        // Add inequality.
        if (addCutFromCoefs(xValue, recourseValue, edgeCoefs, nodeCoefs, RHS,
                            separatedCuts)) {
            nSepCuts++;
        }
    }

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return nSepCuts;
}

void PartialRouteSeparator::setCoefficientsToSet(
    int from, const std::vector<int> &vertices, double coef,
    EdgeValueMap &edgeCoefs) {
    Node u = instance.g.nodeFromId(from);

    for (size_t i = 0; i < vertices.size(); i++) {
        Node v = instance.g.nodeFromId(vertices[i]);
        Edge e = findEdge(instance.g, u, v);
        assert(e != INVALID);
        edgeCoefs[e] = coef;
    }
}

void PartialRouteSeparator::setCoefficientsInsideSet(
    const std::vector<int> &vertices, double coef, EdgeValueMap &edgeCoefs) {
    if (static_cast<int>(vertices.size()) <= 1) {
        return;
    }

    for (size_t i = 0; i < vertices.size(); i++) {
        Node u = instance.g.nodeFromId(vertices[i]);
        for (size_t j = i + 1; j < vertices.size(); j++) {
            Node v = instance.g.nodeFromId(vertices[j]);
            Edge e = findEdge(instance.g, u, v);
            assert(e != INVALID);
            edgeCoefs[e] = coef;
        }
    }
}

void PartialRouteSeparator::setCutCoefficientsNew(
    const PartialRoute &partialRoute, double recourseCost,
    EdgeValueMap &edgeCoefs, double &RHS, bool simpleAdherence, bool improved) {
    // Get number of nodes in the partial route.
    int countPartialRouteNodes = 0;
    int partialRouteSize = static_cast<int>(partialRoute.entries.size());
    for (int i = 0; i < partialRouteSize; i++) {
        countPartialRouteNodes +=
            static_cast<int>(partialRoute.entries[i].vertices.size());
    }
    std::vector<double> alpha(partialRouteSize, 1);
    std::vector<double> beta(partialRouteSize - 1, 1);
    RHS = 2.0 - countPartialRouteNodes;

    // Set alphas and betas.
    if (partialRouteSize == 3) {
        if (!improved ||
            static_cast<int>(partialRoute.entries[0].vertices.size()) == 1 ||
            static_cast<int>(
                partialRoute.entries[partialRouteSize - 1].vertices.size()) ==
                1) {
            alpha[1] += 1;
            RHS +=
                1.0 - static_cast<int>(partialRoute.entries[1].vertices.size());
        }
    } else if (partialRouteSize >= 2) {
        if (!improved ||
            static_cast<int>(partialRoute.entries[0].vertices.size()) == 1) {
            alpha[1] += 1;
            RHS +=
                1.0 - static_cast<int>(partialRoute.entries[1].vertices.size());
        }

        if (!improved ||
            static_cast<int>(
                partialRoute.entries[partialRouteSize - 1].vertices.size()) ==
                1) {
            alpha[partialRouteSize - 2] += 1;
            RHS +=
                1.0 -
                static_cast<int>(
                    partialRoute.entries[partialRouteSize - 2].vertices.size());
        }
    }

    // Set coefficients of arcs incident to the depot.
    if (!simpleAdherence) {
        setCoefficientsToSet(0, partialRoute.entries[0].vertices, -recourseCost,
                             edgeCoefs);
        alpha[0] += 2;
        RHS += -2.0 * static_cast<int>(partialRoute.entries[0].vertices.size());

        if (partialRouteSize >= 2) {
            setCoefficientsToSet(
                0, partialRoute.entries[partialRouteSize - 1].vertices,
                -recourseCost, edgeCoefs);
            alpha[partialRouteSize - 1] += 2;
            beta[0] += 1;
            beta[partialRouteSize - 2] += 1;
            RHS +=
                -2.0 *
                static_cast<int>(
                    partialRoute.entries[partialRouteSize - 1].vertices.size());
        }
    }

    // Set coefficients of the remaining arcs in the partial route.
    for (int i = 0; i < partialRouteSize; i++) {
        if (static_cast<int>(partialRoute.entries[i].vertices.size()) > 1) {
            setCoefficientsInsideSet(partialRoute.entries[i].vertices,
                                     -alpha[i] * recourseCost, edgeCoefs);
        }

        if (i < partialRouteSize - 1) {
            if (static_cast<int>(partialRoute.entries[i].vertices.size()) ==
                1) {
                setCoefficientsToSet(partialRoute.entries[i].vertices[0],
                                     partialRoute.entries[i + 1].vertices,
                                     -beta[i] * recourseCost, edgeCoefs);
            } else {
                setCoefficientsToSet(partialRoute.entries[i + 1].vertices[0],
                                     partialRoute.entries[i].vertices,
                                     -beta[i] * recourseCost, edgeCoefs);
            }
        }
    }

    RHS *= recourseCost;
}

void PartialRouteSeparator::setCutCoefficientsFromDual(
    const PartialRoute &partialRoute, const EdgeValueMap &xValue,
    double recourseCost,
    const std::unordered_map<std::pair<int, int>, double,
                             boost::hash<std::pair<int, int>>> &betaDuals,
    const std::unordered_map<std::tuple<int, int, int>, double,
                             boost::hash<std::tuple<int, int, int>>>
        &alphaDuals,
    EdgeValueMap &edgeCoefs, NodeValueMap &nodeCoefs, double &RHS) const {
    // Zero out coefficients.
    for (const PartialRouteEntry &entry : partialRoute.entries) {
        for (int id : entry.vertices) {
            nodeCoefs[instance.g.nodeFromId(id)] = 0.0;
        }
    }
    RHS = 0.0;

    // Auxilary vector to store accumulated demand.
    std::vector<std::vector<int>> sumDemands(
        instance.nScenarios, std::vector<int>(partialRoute.entries.size(), 0));
    std::vector<int> sumSizes(partialRoute.entries.size(), 0);
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        if (partialRoute.totalScenarioDemands[scenarioId] <=
            instance.capacity) {
            continue;
        }

        int accLoad = 0;
        for (size_t i = 0; i < partialRoute.entries.size(); i++) {
            accLoad += partialRoute.entries[i].scenarioDemands[scenarioId];
            sumDemands[scenarioId][i] = accLoad;
        }
    }

    int accSize = 0;
    for (size_t i = 0; i < partialRoute.entries.size(); i++) {
        accSize += partialRoute.entries[i].vertices.size();
        sumSizes[i] = accSize;
    }

    // Create copies of dual variables to modify later.
    std::vector<std::vector<double>> nodeCoefsAux(
        instance.nScenarios,
        std::vector<double>(instance.n, 0.0)); // [scenarioId][nodeId]
    std::unordered_map<std::pair<int, int>, double,
                       boost::hash<std::pair<int, int>>>
        alphaAux; // [startIdx, endIdx]

    for (const auto &[key, value] : betaDuals) {
        const auto &[scenarioId, id] = key;
        nodeCoefsAux[scenarioId][id] = value;
    }

    // Iterate through alphaDuals to form alphaAux and update betaAux.
    for (const auto &[key, value] : alphaDuals) {
        double alpha = value / static_cast<double>(instance.nScenarios);
        const auto &[scenarioId, startIdx, endIdx] = key;
        assert(partialRoute.totalScenarioDemands[scenarioId] >
               instance.capacity);

        alphaAux[{startIdx, endIdx}] += alpha;
        for (size_t i = startIdx; i < endIdx; i++) {
            for (int id : partialRoute.entries[i].vertices) {
                nodeCoefsAux[scenarioId][id] +=
                    value; // Not divided by nScenarios here.
            }
        }

        // Increment RHS by alpha * (ceil(d^\xi(H') / capacity) - |H'|).
        double startDemand =
            (startIdx == 0) ? 0.0 : sumDemands[scenarioId][startIdx - 1];
        double startSize = (startIdx == 0) ? 0.0 : sumSizes[startIdx - 1];
        RHS += alpha *
               (std::ceil((sumDemands[scenarioId][endIdx - 1] - startDemand) /
                          instance.capacity) -
                (sumSizes[endIdx - 1] - startSize));
    }

    // Set node coefficients.
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        if (partialRoute.totalScenarioDemands[scenarioId] <=
            instance.capacity) {
            continue;
        }

        for (size_t i = 0; i < partialRoute.entries.size(); i++) {
            for (int id : partialRoute.entries[i].vertices) {
                Node v = instance.g.nodeFromId(id);
                std::pair<int, int> key{scenarioId, id};
                double aux =
                    nodeCoefsAux[scenarioId]
                                [id]; // This is confusing, if we replace
                                      // in the formula, we divide and
                                      // multiply by nScenarios.
                nodeCoefs[v] = std::max(
                    nodeCoefs[v], (aux / instance.getEdgeRecourseCost(v)));
                if (betaDuals.find(key) != betaDuals.end()) {
                    RHS += betaDuals.at(key) /
                           static_cast<double>(instance.nScenarios);
                }
            }
        }
    }

    // Set edge coefficients.
    double checkRHS = RHS;
    for (const auto &[key, alpha] : alphaAux) {
        const auto &[startIdx, endIdx] = key;
        std::vector<int> vertices;

        for (size_t i = startIdx; i < endIdx; i++) {
            for (int id : partialRoute.entries[i].vertices) {
                vertices.push_back(id);
            }
        }

        for (size_t i = 0; i < vertices.size(); i++) {
            Node u = instance.g.nodeFromId(vertices[i]);
            for (size_t j = i + 1; j < vertices.size(); j++) {
                Node v = instance.g.nodeFromId(vertices[j]);
                Edge e = findEdge(instance.g, u, v);
                assert(e != INVALID);
                edgeCoefs[e] -= alpha;
                checkRHS += alpha * xValue[e];
            }
        }
    }

    // Check if the cut dominates the ILS cut.
    for (size_t i = 0; i < partialRoute.entries.size(); i++) {
        for (int id : partialRoute.entries[i].vertices) {
            Node v = instance.g.nodeFromId(id);
            assert(nodeCoefs[v] <= 1.0 + EpsForIntegrality);
        }
    }
    assert(checkRHS >= partialRoute.excess * recourseCost - EpsForIntegrality);
}

// The return value of this indicates if the inequality is violated or not.
bool PartialRouteSeparator::addCutFromCoefs(
    const EdgeValueMap &xValue, const NodeValueMap &recourseValue,
    const EdgeValueMap &edgeCoefs, const NodeValueMap &nodeCoefs, double RHS,
    std::vector<CutData> &separatedCuts) {
    separatedCuts.emplace_back(CutData());
    CutData &cutData = separatedCuts.back();
    cutData.RHS = RHS;
    cutData.sense = '>';
    cutData.LHS = 0.0;
    cutData.name = "PartialRouteCut";

    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (std::abs(nodeCoefs[v]) >= EpsForIntegrality) {
            assert(instance.g.id(v) >= 1 && instance.g.id(v) < instance.n);
            cutData.nodePairs.push_back({instance.g.id(v), nodeCoefs[v]});
            cutData.LHS += nodeCoefs[v] * recourseValue[v];
            cutData.customers.push_back(instance.g.id(v));
        }
    }

    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        if (edgeCoefs[e] != 0) {
            cutData.edgePairs.push_back({e, edgeCoefs[e]});
            cutData.LHS += edgeCoefs[e] * xValue[e];
        }
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