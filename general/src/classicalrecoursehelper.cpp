#include "classicalrecoursehelper.h"
#include "svrpinstance.h"

ClassicalRecourseHelper::ClassicalRecourseHelper(const SVRPInstance &instance_)
    : instance(instance_) {}

ClassicalRecourseHelper::~ClassicalRecourseHelper() {}

double ClassicalRecourseHelper::getRecourseCost(
    const std::vector<std::vector<int>> &routes) const {
    double totalCost = 0.0;

    for (size_t i = 0; i < routes.size(); i++) {
        totalCost += getRouteRecourseCost(routes[i]);
    }

    return totalCost;
}

double ClassicalRecourseHelper::getRouteRecourseCost(
    const std::vector<int> &route) const {
    double totalCost1 = 0.0;
    double totalCost2 = 0.0;

    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        totalCost1 += getRouteRecourseCostInScenario(route, scenarioId, false);
        totalCost2 += getRouteRecourseCostInScenario(route, scenarioId, true);
    }

    return std::min(totalCost1, totalCost2) /
           static_cast<double>(instance.nScenarios);
}

double ClassicalRecourseHelper::getRouteRecourseCostInScenario(
    const std::vector<int> &route, int scenarioId, bool reverse) const {
    if (route.empty()) {
        return 0.0;
    }

    double cost = 0.0;
    double curr = 0.0;
    int start = (reverse ? route.size() - 1 : 0);
    int end = (reverse ? -1 : route.size());

    for (int i = start; i != end; reverse ? --i : ++i) {
        assert(route[i] >= 0 && route[i] < instance.n);
        curr += instance.scenariosMatrix[route[i]][scenarioId];

        if (curr > instance.capacity) {
            cost +=
                instance.getEdgeRecourseCost(instance.g.nodeFromId(route[i]));
            curr -= instance.capacity;
        }
    }

    return cost;
}

double ClassicalRecourseHelper::getRecourseCostInScenario(
    const std::vector<std::vector<int>> &routes, int scenarioId) const {
    assert(scenarioId < instance.nScenarios);
    double totalCost = 0.0;

    for (auto &route : routes) {
        double cost1 = getRouteRecourseCostInScenario(route, scenarioId, false);
        double cost2 = getRouteRecourseCostInScenario(route, scenarioId, true);
        totalCost += std::min(cost1, cost2);
    }

    return totalCost;
}

// Check if the route is in the wrong direction.
bool ClassicalRecourseHelper::isRouteReversed(
    const std::vector<int> &route) const {
    double totalCost1 = 0.0;
    double totalCost2 = 0.0;

    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        double cost1 = 0.0;
        double cost2 = 0.0;
        double curr1 = 0.0;
        double curr2 = 0.0;
        size_t i, j;
        Edge e1, e2;

        for (i = 0; i < route.size(); i++) {
            j = route.size() - i - 1;

            assert(route[i] >= 1 && route[i] < instance.n);
            assert(route[j] >= 1 && route[j] < instance.n);
            curr1 += instance.scenariosMatrix[route[i]][scenarioId];
            curr2 += instance.scenariosMatrix[route[j]][scenarioId];

            if (curr1 > instance.capacity) {
                cost1 += instance.getEdgeRecourseCost(
                    instance.g.nodeFromId(route[i]));
                curr1 -= instance.capacity;
            }

            if (curr2 > instance.capacity) {
                cost2 += instance.getEdgeRecourseCost(
                    instance.g.nodeFromId(route[j]));
                curr2 -= instance.capacity;
            }
        }

        totalCost1 += cost1;
        totalCost2 += cost2;
    }

    return totalCost2 < totalCost1;
}

double ClassicalRecourseHelper::getPartialRouteRecourseAndFailureNodes(
    const PartialRoute &partialRoute, NodeValueMap &nodeRecourseMap) const {
    NodeValueMap nodeRecourseMap1(instance.g, 0.0);
    NodeValueMap nodeRecourseMap2(instance.g, 0.0);
    double recourseCost1 = 0.0;
    double recourseCost2 = 0.0;
    int routeSize = static_cast<int>(partialRoute.entries.size());

    for (int s = 0; s < instance.nScenarios; s++) {
        int load1 = 0.0;
        int load2 = 0.0;

        for (int i = 0; i < routeSize; i++) {
            // Compute the partial route entry demand.
            int demand1 = 0;
            int demand2 = 0;
            for (int id : partialRoute.entries[i].vertices) {
                demand1 += instance.scenariosMatrix[id][s];
            }
            for (int id : partialRoute.entries[routeSize - 1 - i].vertices) {
                demand2 += instance.scenariosMatrix[id][s];
            }

            // Get fail count.
            int fail1 = 0;
            int fail2 = 0;
            if (load1 == 0) {
                fail1 = static_cast<int>(
                    std::ceil(demand1 / instance.capacity) - 1);
            } else {
                int remainder = load1 % static_cast<int>(instance.capacity);
                fail1 = static_cast<int>(
                    std::ceil((demand1 + remainder) / instance.capacity) -
                    std::ceil(remainder / instance.capacity));
            }

            if (load2 == 0) {
                fail2 = static_cast<int>(
                    std::ceil(demand2 / instance.capacity) - 1);
            } else {
                int remainder = load2 % static_cast<int>(instance.capacity);
                fail2 = static_cast<int>(
                    std::ceil((demand2 + remainder) / instance.capacity) -
                    std::ceil(remainder / instance.capacity));
            }

            // Update recourse cost.
            if (fail1 >= 1) {
                assert(fail1 <= partialRoute.entries[i].vertices.size());
                for (int j = 0; j < fail1; j++) {
                    Node v = instance.g.nodeFromId(
                        partialRoute.entries[i].vertices[j]);
                    recourseCost1 += instance.getEdgeRecourseCost(v);
                }

                if (partialRoute.isRoute) {
                    assert(partialRoute.entries[i].vertices.size() == 1);
                    int id = partialRoute.entries[i].vertices[0];
                    Node v = instance.g.nodeFromId(id);
                    nodeRecourseMap1[v] +=
                        instance.getEdgeRecourseCost(v) /
                        static_cast<double>(instance.nScenarios);
                }
            }

            if (fail2 >= 1) {
                assert(fail2 <=
                       partialRoute.entries[routeSize - 1 - i].vertices.size());
                for (int j = 0; j < fail2; j++) {
                    Node v = instance.g.nodeFromId(
                        partialRoute.entries[routeSize - 1 - i].vertices[j]);
                    recourseCost2 += instance.getEdgeRecourseCost(v);
                }

                if (partialRoute.isRoute) {
                    assert(partialRoute.entries[routeSize - 1 - i]
                               .vertices.size() == 1);
                    int id =
                        partialRoute.entries[routeSize - 1 - i].vertices[0];
                    Node v = instance.g.nodeFromId(id);
                    nodeRecourseMap2[v] +=
                        instance.getEdgeRecourseCost(v) /
                        static_cast<double>(instance.nScenarios);
                }
            }

            load1 += demand1;
            load2 += demand2;
        }
    }

    // Set count failures.
    if (partialRoute.isRoute) {
        if (recourseCost1 <= recourseCost2) {
            for (NodeIt v(instance.g); v != INVALID; ++v) {
                nodeRecourseMap[v] = nodeRecourseMap1[v];
            }
        } else {
            for (NodeIt v(instance.g); v != INVALID; ++v) {
                nodeRecourseMap[v] = nodeRecourseMap2[v];
            }
        }
    }

    return std::min(recourseCost1, recourseCost2) /
           static_cast<double>(instance.nScenarios);
}

double ClassicalRecourseHelper::getRouteRecourseCostAndFailures(
    const std::vector<int> &route, NodeVectorValueMap &failures) const {
    if (route.empty()) {
        return 0.0;
    }

    NodeVectorValueMap failures1(instance.g,
                                 std::vector<double>(instance.nScenarios, 0.0));
    NodeVectorValueMap failures2(instance.g,
                                 std::vector<double>(instance.nScenarios, 0.0));
    double cost1 = 0.0;
    double cost2 = 0.0;
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        double curr1 = 0.0;
        double curr2 = 0.0;

        for (size_t i = 0; i < route.size(); i++) {
            size_t j = route.size() - 1 - i;
            int id1 = route[i];
            int id2 = route[j];
            assert(id1 >= 1 && id1 < instance.n);
            assert(id2 >= 1 && id2 < instance.n);
            curr1 += instance.scenariosMatrix[id1][scenarioId];
            curr2 += instance.scenariosMatrix[id2][scenarioId];

            if (curr1 > instance.capacity) {
                Node v = instance.g.nodeFromId(id1);
                cost1 += instance.getEdgeRecourseCost(v);
                curr1 -= instance.capacity;
                failures1[v][scenarioId] += 1;
            }

            if (curr2 > instance.capacity) {
                Node v = instance.g.nodeFromId(id2);
                cost2 += instance.getEdgeRecourseCost(v);
                curr2 -= instance.capacity;
                failures2[v][scenarioId] += 1;
            }
        }
    }

    NodeVectorValueMap &chosen = (cost1 <= cost2 ? failures1 : failures2);
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        for (int scenarioId = 0; scenarioId < instance.nScenarios;
             scenarioId++) {
            failures[v][scenarioId] = chosen[v][scenarioId];
        }
    }

    return std::min(cost1, cost2) / static_cast<double>(instance.nScenarios);
}