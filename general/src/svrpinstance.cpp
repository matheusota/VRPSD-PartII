#include "svrpinstance.h"
#include "optimalrecoursehelper.h"

SVRPInstance::SVRPInstance()
    : posx(g), posy(g), demand(g, 0.0), weight(g, 0.0), upperBoundMap(g, 1.0),
      d_posx(d_g), d_posy(d_g), d_demand(d_g, 0.0), d_weight(d_g, 0.0),
      classicalRecourseHelper(*this), optimalRecourseHelper(*this) {
    n = m = k = depot = nrows = 0;
    capacity = 0.0;
    nScenarios = 0;
}

void SVRPInstance::sortScenarios() {
    for (int scenarioId = 0; scenarioId < nScenarios; scenarioId++) {
        double totalDemand = 0.0;

        for (int i = 1; i < n; i++) {
            totalDemand += scenariosMatrix[i][scenarioId];
        }

        sortedScenarios.push_back(std::make_pair(-totalDemand, scenarioId));
    }

    std::sort(sortedScenarios.begin(), sortedScenarios.end());
}

Node SVRPInstance::getNodeFromDNode(DNode v) const {
    return g.nodeFromId(d_g.id(v));
}
DNode SVRPInstance::getDNodeFromNode(Node v) const {
    return d_g.nodeFromId(g.id(v));
}

Edge SVRPInstance::getEdgeFromArc(Arc a) const {
    return findEdge(g, getNodeFromDNode(d_g.source(a)),
                    getNodeFromDNode(d_g.target(a)));
}

std::pair<Arc, Arc> SVRPInstance::getArcsFromEdge(Edge e) const {
    return std::make_pair(
        findArc(d_g, getDNodeFromNode(g.u(e)), getDNodeFromNode(g.v(e))),
        findArc(d_g, getDNodeFromNode(g.v(e)), getDNodeFromNode(g.u(e))));
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

double SVRPInstance::getArcRecourseCost(DNode v) const {
    assert(v != INVALID);

    if (d_g.id(v) != depot) {
        Arc depotArc = findArc(d_g, d_g.nodeFromId(depot), v);
        assert(depotArc != INVALID);
        return 2.0 * d_weight[depotArc];
    } else {
        return 0.0;
    }
}