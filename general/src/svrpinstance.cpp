#include "svrpinstance.h"
#include "optimalrecoursehelper.h"

SVRPInstance::SVRPInstance()
    : posx(g), posy(g), demand(g, 0.0), weight(g, 0.0),
      classicalRecourseHelper(*this), optimalRecourseHelper(*this) {
    n = m = k = depot = 0;
    capacity = 0.0;
    nScenarios = 0;
}

void SVRPInstance::print() const {
    std::cout << "Input Description:" << endl;
    std::cout << "n          : " << n << endl;
    std::cout << "m          : " << m << endl;
    std::cout << "k          : " << k << endl;
    std::cout << "capacity   : " << capacity << endl;
    std::cout << "depot      : " << depot << endl;

    std::cout << "nodes     :" << endl;
    for (NodeIt v(g); v != INVALID; ++v) {
        std::cout << " " << g.id(v) << ": " << demand[v] << endl;
    }
    std::cout << endl;

    std::cout << "-------------------------" << endl;
    for (int j = 0; j < nScenarios; j++) {
        std::cout << "SCENARIO " << j << " :";
        for (int i = 0; i < n; i++) {
            std::cout << scenariosMatrix[i][j] << " ";
        }
        std::cout << endl;
    }
}

void SVRPInstance::sortScenarios() {
    for (int scenarioId = 0; scenarioId < nScenarios; scenarioId++) {
        double totalDemand = 0.0;

        for (int i = 1; i < n; i++) {
            totalDemand += scenariosMatrix[i][scenarioId];
        }

        sortedScenarios.push_back({-totalDemand, scenarioId});
    }

    std::sort(sortedScenarios.begin(), sortedScenarios.end());
}

// Check the feasibility according to the average demands.
bool SVRPInstance::checkDeterministicRoute(
    const std::vector<int> &route) const {
    double load = 0.0;
    for (int i : route) {
        load += demand[g.nodeFromId(i)];
        if (load > capacity) {
            return false;
        }
    }

    return true;
}

double SVRPInstance::getEdgeRecourseCost(Node v) const {
    assert(v != INVALID);

    if (g.id(v) != depot) {
        Edge depotEdge = findEdge(g, g.nodeFromId(depot), v);
        assert(depotEdge != INVALID);
        return 2.0 * weight[depotEdge];
    } else {
        return 0.0;
    }
}