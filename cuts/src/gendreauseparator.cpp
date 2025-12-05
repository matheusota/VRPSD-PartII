#include "gendreauseparator.h"

GendreauSeparator::GendreauSeparator(const SVRPInstance &instance,
                                     const Params &params)
    : instance(instance), params(params) {}

int GendreauSeparator::separateGendreauCutsPerScenario(
    const EdgeValueMap &xValue, const std::vector<double> &scenariosRecourse,
    std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    std::vector<std::vector<int>> routes;

    // Only separate Gendreau cuts for integer solutions.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        if (std::abs(xValue[e] - std::round(xValue[e])) >= EpsForIntegrality) {
            return 0;
        }
    }

    if (!captureRoutes(xValue, routes)) {
        std::cout << "ERROR: GendreauSeparator cannot capture routes."
                  << std::endl;
        exit(1);
    }

    // For each scenario, compute the recourse.
    int sepCuts = 0;
    for (size_t scenarioId = 0; scenarioId < scenariosRecourse.size();
         scenarioId++) {
        sepCuts += addGendreauCutForScenario(xValue, scenariosRecourse, routes,
                                             scenarioId, separatedCuts);
    }

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return sepCuts;
}

int GendreauSeparator::separateAggregatedGendreauCut(
    const EdgeValueMap &xValue, double totalRecourse,
    std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    // Only separate Gendreau cuts for integer solutions.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        if (std::abs(xValue[e] - std::round(xValue[e])) >= EpsForIntegrality) {
            return 0;
        }
    }

    std::vector<std::vector<int>> routes;
    if (!captureRoutes(xValue, routes)) {
        std::cout << "ERROR: GendreauSeparator cannot capture routes."
                  << std::endl;
        exit(1);
    }

    int sepCuts =
        addAggregatedGendreauCut(xValue, totalRecourse, routes, separatedCuts);
    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return sepCuts;
}

// Add a cut of the form y >= Q(x) (\sum_{e \in \bar{E}} x_e - n + k + 2),
// as described in Gendreau (1995). This add a single cut for each scenario.
int GendreauSeparator::addGendreauCutForScenario(
    const EdgeValueMap &xValue, const std::vector<double> &scenariosRecourse,
    const std::vector<std::vector<int>> &routes, int scenarioId,
    std::vector<CutData> &separatedCuts) {
    double currRecourseCost = scenariosRecourse[scenarioId];
    double recourseCost =
        instance.classicalRecourseHelper.getRecourseCostInScenario(routes,
                                                                   scenarioId);
    if (recourseCost <= currRecourseCost + EpsForIntegrality) {
        return 0;
    }

    separatedCuts.emplace_back(CutData());
    CutData &cutData = separatedCuts.back();
    cutData.RHS = recourseCost * (2 + instance.k - instance.n);
    cutData.sense = '>';
    cutData.name = "GendreauCut";
    cutData.LHS = currRecourseCost;
    cutData.scenarioPairs.push_back(std::make_pair(scenarioId, 1.0));

    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) != instance.depot) {
            cutData.nodePairs.push_back({instance.g.id(v), 1.0});
        }
    }

    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        if (xValue[e] >= 1 - EpsForIntegrality &&
            instance.g.id(instance.g.u(e)) != instance.depot &&
            instance.g.id(instance.g.v(e)) != instance.depot) {
            cutData.edgePairs.push_back({e, -recourseCost});
            cutData.LHS -= recourseCost * xValue[e];
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

// Add a cut of the form y >= Q(x) (\sum_{e \in \bar{E}} x_e - n + k + 2),
// as described in Gendreau (1995). This aggregates the scenarios.
int GendreauSeparator::addAggregatedGendreauCut(
    const EdgeValueMap &xValue, double totalRecourse,
    const std::vector<std::vector<int>> &routes,
    std::vector<CutData> &separatedCuts) {

    double currRecourse =
        instance.classicalRecourseHelper.getRecourseCost(routes);

    if (currRecourse <= totalRecourse + EpsForIntegrality) {
        return 0;
    }

    separatedCuts.emplace_back(CutData());
    CutData &cutData = separatedCuts.back();
    cutData.RHS = currRecourse * (2 + instance.k - instance.n);
    cutData.sense = '>';
    cutData.name = "GendreauCut";
    cutData.LHS = totalRecourse;

    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        if (xValue[e] >= 1 - EpsForIntegrality &&
            instance.g.id(instance.g.u(e)) != instance.depot &&
            instance.g.id(instance.g.v(e)) != instance.depot) {
            cutData.edgePairs.push_back(std::make_pair(e, -currRecourse));
            cutData.LHS -= currRecourse * xValue[e];
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

// This checks if an integer solution is valid for the CVRP and capture each of
// the routes in the solution.
bool GendreauSeparator::captureRoutes(const EdgeValueMap &xValue,
                                      std::vector<std::vector<int>> &routes) {
    std::vector<int> currRoute;
    Node depot = instance.g.nodeFromId(instance.depot);
    Node curr = depot;
    Node next;
    bool hasNext = true;
    int count = 1;
    int load = 0.0;
    EdgeValueMap edgeCount(instance.g);
    lemon::mapCopy(instance.g, xValue, edgeCount);

    // Capture the routes.
    while (hasNext) {
        hasNext = false;

        // get next node
        for (IncEdgeIt e(instance.g, curr); e != INVALID; ++e) {
            if (edgeCount[e] > EpsForIntegrality) {
                next = (instance.g.u(e) == curr) ? instance.g.v(e)
                                                 : instance.g.u(e);
                edgeCount[e] -= 1.0;
                hasNext = true;
                break;
            }
        }

        // go to next vertex and update.
        if (hasNext) {
            // back to depot
            if (next == depot) {
                routes.push_back(currRoute);
                currRoute.clear();
                load = 0.0;
            }
            // adding vertex to curr route
            else {
                load += instance.demand[next];
                count++;
                currRoute.push_back(instance.g.id(next));

                // capacity exceeded
                if (load > instance.capacity) {
                    return false;
                }
            }

            curr = next;
        }
    }

    // not all vertices were visited, solution is not feasible
    if (count < instance.n) {
        return false;
    } else {
        return true;
    }
}