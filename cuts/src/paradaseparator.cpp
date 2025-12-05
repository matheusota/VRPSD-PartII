#include "paradaseparator.h"

ParadaSeparator::ParadaSeparator(const SVRPInstance &instance,
                                 const Params &params)
    : instance(instance), params(params) {}

int ParadaSeparator::separateSetCut(const EdgeValueMap &xValue,
                                    const NodeValueMap &recourseValue,
                                    const std::vector<int> &customers,
                                    std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    int sepCuts = addSetCut(xValue, recourseValue, customers, separatedCuts);

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return sepCuts;
}

void ParadaSeparator::buildSupportGraph(const EdgeValueMap &xValue,
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

int ParadaSeparator::separateParadaCuts(const EdgeValueMap &xValue,
                                        const NodeValueMap &recourseValue,
                                        std::vector<CutData> &separatedCuts,
                                        bool isFeasible) {
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

        int setCuts = 0;
        if (params.paradaSet) {
            setCuts = addSetCut(xValue, recourseValue, componentCustomers[i],
                                separatedCuts);
        }

        // Component is a path.
        int pathCuts = 0;
        if (componentCustomers[i].size() > 1 && setCuts == 0 &&
            componentNonZero[i] ==
                static_cast<int>(componentCustomers[i].size()) - 1) {
            pathCuts = collectPathCuts(supportGraph, xValue, recourseValue,
                                       separatedCuts, componentCustomers[i],
                                       isFeasible);
        }
        nSepCuts += setCuts + pathCuts;
    }

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return nSepCuts;
}

int ParadaSeparator::collectPathCuts(const Graph &supportGraph,
                                     const EdgeValueMap &xValue,
                                     const NodeValueMap &recourseValue,
                                     std::vector<CutData> &separatedCuts,
                                     const std::vector<int> &customers,
                                     bool isFeasible) {
    // Compute the degree of each customer.
    Node first = supportGraph.nodeFromId(customers[0]);

    // Search for a customer with degree one (assume we have more than 1
    // customer).
    for (size_t i = 0; i < customers.size(); i++) {
        Node v = supportGraph.nodeFromId(customers[i]);
        int degree = 0;
        for (IncEdgeIt e(supportGraph, v); e != INVALID; ++e) {
            degree++;
        }

        // Component is not a path!
        if (degree >= 3) {
            return 0;
        }

        if (degree == 1) {
            first = v;
            break;
        }
    }

    // Starting from a customer with degree 1, construct the route.
    Node prev = INVALID;
    Node curr = first;
    bool hasNext = true;
    std::vector<int> route;
    while (hasNext) {
        hasNext = false;
        route.push_back(supportGraph.id(curr));

        for (IncEdgeIt e(supportGraph, curr); e != INVALID; ++e) {
            Node v = (supportGraph.u(e) == curr) ? supportGraph.v(e)
                                                 : supportGraph.u(e);
            if (v != prev) {
                prev = curr;
                curr = v;
                hasNext = true;
                break;
            }
        }
    }

    // Add path cuts and generate subroutes removing at most 5 customers.
    int addedCuts = 0;
    addedCuts += addPathCut(route, xValue, recourseValue, separatedCuts);

    // if (isFeasible) {
    //     int nDelete = std::min(static_cast<int>(route.size()) - 1, 5);

    //     for (int j = 1; j <= nDelete; j++) {
    //         // Create a vector with last j entries set to false.
    //         std::vector<bool> permutation(route.size(), false);
    //         for (size_t t = route.size() - j; t < route.size(); t++) {
    //             permutation[t] = true;
    //         }

    //         // Generate permutations to get all subroutes.
    //         do {
    //             std::vector<int> newRoute;
    //             for (size_t i = 0; i < route.size(); i++) {
    //                 if (!permutation[i]) {
    //                     newRoute.push_back(route[i]);
    //                 }
    //             }

    //             addedCuts +=
    //                 addPathCut(newRoute, xValue, recourseValue,
    //                 separatedCuts);
    //         } while (
    //             std::next_permutation(permutation.begin(),
    //             permutation.end()));
    //     }
    // }

    return addedCuts;
}

// The return value of this indicates if the inequality is violated or not.
int ParadaSeparator::addPathCut(const std::vector<int> &route,
                                const EdgeValueMap &xValue,
                                const NodeValueMap &recourseValue,
                                std::vector<CutData> &separatedCuts) {
    // Check route edges first.
    double edgeSumCheck = 0.0;
    for (size_t j = 0; j < route.size() - 1; j++) {
        Edge e = findEdge(instance.g, instance.g.nodeFromId(route[j]),
                          instance.g.nodeFromId(route[j + 1]));
        assert(e != INVALID);
        edgeSumCheck += 1.0 - xValue[e];
        if (edgeSumCheck >= 1.0) {
            return 0;
        }
    }

    // Compute recourse.
    double recourseCost =
        instance.optimalRecourseHelper.getRouteRecourseCost(route);
    if (recourseCost <= EpsForIntegrality) {
        return 0;
    }

    // Construct cut.
    separatedCuts.emplace_back(CutData());
    CutData &cutData = separatedCuts.back();
    cutData.RHS = recourseCost * (2.0 - static_cast<double>(route.size()));
    cutData.sense = '>';
    cutData.LHS = 0.0;
    cutData.name = "PathCut";

    for (int id : route) {
        Node v = instance.g.nodeFromId(id);
        assert(id != instance.depot);
        cutData.nodePairs.push_back(std::make_pair(id, 1.0));
        cutData.LHS += recourseValue[v];
        cutData.customers.push_back(id);
    }

    for (size_t j = 0; j < route.size() - 1; j++) {
        Edge e = findEdge(instance.g, instance.g.nodeFromId(route[j]),
                          instance.g.nodeFromId(route[j + 1]));
        assert(e != INVALID);
        cutData.edgePairs.push_back(std::make_pair(e, -recourseCost));
        cutData.LHS -= recourseCost * xValue[e];
    }

    if (cutData.LHS >= cutData.RHS - EpsForIntegrality) {
        separatedCuts.pop_back();
        return 0;
    } else {
        cutData.violation = cutData.RHS - cutData.LHS;
        nPathCuts++;
        return 1;
    }
}

// The return value of this indicates if the inequality is violated or not.
int ParadaSeparator::addSetCut(const EdgeValueMap &xValue,
                               const NodeValueMap &recourseValue,
                               const std::vector<int> &customers,
                               std::vector<CutData> &separatedCuts) {
    double totalDemand = 0.0;
    for (int id : customers) {
        Node v = instance.g.nodeFromId(id);
        assert(id >= 1 && id < instance.n && v != INVALID);
        totalDemand += instance.demand[v];
    }
    int nVehicles =
        (params.basicVRPSD)
            ? 1
            : static_cast<int>(std::ceil(totalDemand / instance.capacity));

    // Compute recourse.
    double lowerBound = instance.optimalRecourseHelper.getRecourseLowerBound(
        customers, nVehicles);
    if (lowerBound <= EpsForIntegrality) {
        return 0;
    }

    // Construct cut.
    separatedCuts.emplace_back(CutData());
    CutData &cutData = separatedCuts.back();
    cutData.name = "SetCut";
    cutData.RHS =
        lowerBound * (nVehicles + 1 - static_cast<int>(customers.size()));
    cutData.sense = '>';
    cutData.LHS = 0.0;

    for (int id : customers) {
        Node v = instance.g.nodeFromId(id);
        assert(v != INVALID);
        cutData.nodePairs.push_back(std::make_pair(instance.g.id(v), 1.0));
        cutData.LHS += recourseValue[v];
        cutData.customers.push_back(id);
    }

    for (size_t j = 0; j < customers.size(); j++) {
        for (size_t k = j + 1; k < customers.size(); k++) {
            Edge e = findEdge(instance.g, instance.g.nodeFromId(customers[j]),
                              instance.g.nodeFromId(customers[k]));
            assert(e != INVALID);
            cutData.edgePairs.push_back(std::make_pair(e, -lowerBound));
            cutData.LHS -= lowerBound * xValue[e];
        }
    }

    if (cutData.LHS >= cutData.RHS - EpsForIntegrality) {
        separatedCuts.pop_back();
        return 0;
    } else {
        cutData.violation = cutData.RHS - cutData.LHS;
        nSetCuts++;
        return 1;
    }
}