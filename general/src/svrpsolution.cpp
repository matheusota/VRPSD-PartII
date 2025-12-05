#include "svrpsolution.h"

void SVRPSolution::setRoutesFromSolution(const SVRPInstance &instance,
                                         EdgeIntMap &edgeCount) {
    std::vector<int> currRoute;
    Node depot = instance.g.nodeFromId(0);
    Node curr = depot;
    Node next;
    bool hasNext = true;
    routes.clear();

    while (hasNext) {
        hasNext = false;

        // get next node
        for (IncEdgeIt e(instance.g, curr); e != INVALID; ++e) {
            if (edgeCount[e] > 0) {
                next = (instance.g.u(e) == curr) ? instance.g.v(e)
                                                 : instance.g.u(e);
                edgeCount[e]--;
                hasNext = true;
                break;
            }
        }

        if (hasNext) {
            // back to depot
            if (next == depot) {
                routes.push_back(currRoute);
                currRoute.clear();
            }
            // adding vertex to curr route
            else {
                currRoute.push_back(instance.g.id(next));
            }

            curr = next;
        }
    }

    // fix route orientations
    for (auto route : routes) {
        if (instance.classicalRecourseHelper.isRouteReversed(route)) {
            std::reverse(route.begin(), route.end());
        }
    }
}

void SVRPSolution::setRoutesFromSolution(const SVRPInstance &instance,
                                         ArcIntMap &edgeCount) {
    std::vector<int> currRoute;
    DNode depot = instance.d_g.nodeFromId(0);
    DNode curr = depot;
    DNode next;
    bool hasNext = true;
    routes.clear();

    while (hasNext) {
        hasNext = false;

        // get next node
        for (OutArcIt e(instance.d_g, curr); e != INVALID; ++e) {
            if (edgeCount[e] > 0) {
                next = instance.d_g.target(e);
                edgeCount[e]--;
                hasNext = true;
                break;
            }
        }

        if (hasNext) {
            // back to depot
            if (next == depot) {
                routes.push_back(currRoute);
                currRoute.clear();
            }
            // adding vertex to curr route
            else {
                currRoute.push_back(instance.d_g.id(next));
            }

            curr = next;
        }
    }
}

std::vector<Edge>
SVRPSolution::getEdgesInSolution(const SVRPInstance &instance) const {
    std::vector<Edge> edgesInRoutes;

    for (auto &route : routes) {
        int prev = instance.depot;

        for (int curr : route) {
            edgesInRoutes.push_back(findEdge(instance.g,
                                             instance.g.nodeFromId(prev),
                                             instance.g.nodeFromId(curr)));

            prev = curr;
        }

        edgesInRoutes.push_back(
            findEdge(instance.g, instance.g.nodeFromId(prev),
                     instance.g.nodeFromId(instance.depot)));
    }

    return edgesInRoutes;
}