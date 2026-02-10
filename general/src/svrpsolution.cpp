#include "svrpsolution.h"

void SVRPSolution::print() const {
    std::cout << "Solved!" << std::endl;
    std::cout << "Solution Cost: " << cost << std::endl;
    std::cout << "Optimal Recourse Cost: " << optimalRecourseCost << std::endl;
    std::cout << "Recourse Cost: " << recourseCost << std::endl;
    std::cout << "Time: " << time << "s" << std::endl;
    std::cout << "Nodes explored: " << nodesExplored << std::endl;
    std::cout << "CVRPSEP Time: " << cvrpsepTime << std::endl;
    std::cout << "Parada Time: " << paradaTime << std::endl;
    std::cout << "Partial Route Time: " << partialRouteTime << std::endl;
    std::cout << "SRI Time: " << sriTime << std::endl;
    std::cout << "Aggregated SRI Time: " << aggregatedSriTime << std::endl;
    std::cout << "Routes:" << std::endl;
    for (const auto &route : routes) {
        for (int i : route) {
            std::cout << std::to_string(i) << " ";
        }
        std::cout << std::endl;
    }
}

void SVRPSolution::save(std::string outputName) const {
    std::ofstream myfile;
    myfile.open(outputName);

    if (!myfile) {
        cerr << "Error writing solution to " << outputName << endl;
        exit(0);
    }

    if (upperBound >= 1e6) {
        myfile << "No solution found." << endl;
        return;
    }

    for (auto route : routes) {
        for (int i : route) {
            myfile << to_string(i) << " ";
        }
        myfile << endl;
    }

    myfile.close();
}

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

        // Get the next node.
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
            // Back to the depot.
            if (next == depot) {
                routes.push_back(currRoute);
                currRoute.clear();
            }
            // Add vertex to curr route.
            else {
                currRoute.push_back(instance.g.id(next));
            }

            curr = next;
        }
    }

    // Fix route orientations.
    for (auto route : routes) {
        if (instance.classicalRecourseHelper.isRouteReversed(route)) {
            std::reverse(route.begin(), route.end());
        }
    }

    // Check that all edges were used.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        assert(edgeCount[e] == 0);
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

void SVRPSolution::checkCVRPFeasible(const SVRPInstance &instance) {
    std::vector<bool> visited(instance.n, false);

    for (auto route : routes) {
        double load = 0.0;

        for (int id : route) {
            Node v = instance.g.nodeFromId(id);
            load += instance.demand[v];
            if (visited[id]) {
                cvrpFeasible = false;
                return;
            }
            visited[id] = true;
        }

        if (load >= instance.capacity + 1e-9) {
            cvrpFeasible = false;
            return;
        }
    }

    visited[0] = true;
    for (int i = 0; i < instance.n; i++) {
        if (!visited[i]) {
            cvrpFeasible = false;
            return;
        }
    }

    cvrpFeasible = true;
}