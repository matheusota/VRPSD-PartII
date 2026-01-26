#include "svrpinstance.h"
#include "partialroute.h"

PartialRoute::PartialRoute(const SVRPInstance &instance) : instance(instance) {
    totalScenarioDemands.reserve(instance.nScenarios);
}

void PartialRoute::buildFromRoute(const std::vector<int> &route) {
    for (const int &id : route) {
        Node v = instance.g.nodeFromId(id);
        assert(v != INVALID);
        entries.push_back(
            PartialRouteEntry({id}, instance.getEdgeRecourseCost(v), 0.0, 0.0,
                              0.0, instance.nScenarios));
    }

    // Fill in the scenario demands.
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        double totalScenarioDemand = 0.0;

        for (size_t i = 0; i < entries.size(); i++) {
            int id = entries[i].vertices[0];
            entries[i].scenarioDemands.push_back(
                instance.scenariosMatrix[id][scenarioId]);
            totalScenarioDemand += instance.scenariosMatrix[id][scenarioId];
        }

        totalScenarioDemands.push_back(totalScenarioDemand);
    }
}