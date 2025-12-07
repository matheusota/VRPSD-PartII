#include "partialrouteseparator.h"
#include "partialroutesbuilder.h"

PartialRouteSeparator::PartialRouteSeparator(const SVRPInstance &instance,
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
        // Get partial route costs.
        NodeValueMap nodeRecourseMap(instance.g, 0.0);
        double partialRouteClassicalCost =
            instance.classicalRecourseHelper
                .getPartialRouteRecourseAndFailureNodes(partialRoute,
                                                        nodeRecourseMap);
        double partialRouteOptimalCost = 0.0;
        if (params.scenarioOptimalPRCuts) {
            partialRouteOptimalCost =
                instance.optimalRecourseHelper.getPartialRouteRecourseCost(
                    partialRoute);
        }
        double currCost = 0.0;
        std::vector<Node> nodesToConsider;

        // Pass through all nodes in the partial route.
        for (const PartialRouteEntry &entry : partialRoute.entries) {
            for (int id : entry.vertices) {
                Node v = instance.g.nodeFromId(id);
                currCost += recourseValue[v];
                nodesToConsider.push_back(v);
            }
        }

        // Compare costs.
        bool useSimpleAdherenceCuts = false;
        double recourseCostToUse = partialRouteClassicalCost;
        if (currCost + EpsForIntegrality <= partialRouteOptimalCost) {
            recourseCostToUse = partialRouteOptimalCost;
            useSimpleAdherenceCuts = true;
        }

        if (currCost >= recourseCostToUse - EpsForIntegrality) {
            continue;
        }

        // Setup the coefficients for the inequality.
        double RHS = 0.0;
        EdgeValueMap edgeCoefs(instance.g, 0.0);
        setCutCoefficientsNew(partialRoute, recourseCostToUse, edgeCoefs, RHS,
                              useSimpleAdherenceCuts,
                              params.scenarioOptimalPRCuts);

        // Add inequality.
        if (addCutFromCoefs(xValue, nodesToConsider, recourseValue, edgeCoefs,
                            RHS, separatedCuts)) {
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

// The return value of this indicates if the inequality is violated or not.
bool PartialRouteSeparator::addCutFromCoefs(
    const EdgeValueMap &xValue, const std::vector<Node> &nodesToConsider,
    const NodeValueMap &recourseValue, const EdgeValueMap &edgeCoefs,
    double RHS, std::vector<CutData> &separatedCuts) {
    separatedCuts.emplace_back(CutData());
    CutData &cutData = separatedCuts.back();
    cutData.RHS = RHS;
    cutData.sense = '>';
    cutData.LHS = 0.0;
    cutData.name = "PartialRouteCut";

    for (Node v : nodesToConsider) {
        assert(instance.g.id(v) >= 1 && instance.g.id(v) < instance.n);
        cutData.nodePairs.push_back({instance.g.id(v), 1.0});
        cutData.LHS += recourseValue[v];
        cutData.customers.push_back(instance.g.id(v));
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