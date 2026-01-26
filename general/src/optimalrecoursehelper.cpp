#include "optimalrecoursehelper.h"
#include "svrpinstance.h"
#include "partialroute.h"

OptimalRecourseHelper::OptimalRecourseHelper(const SVRPInstance &instance_)
    : instance(instance_) {
    env.set(GRB_IntParam_LogToConsole, 0);
}

OptimalRecourseHelper::~OptimalRecourseHelper() {}

// Lower bound if `customers` is attended using `nVehicles` vehicles.
double
OptimalRecourseHelper::getRecourseLowerBound(const std::vector<int> &customers,
                                             int nVehicles) const {
    if (customers.empty()) {
        return 0.0;
    }

    // Sort customers by their recourse costs.
    std::vector<double> sortedRecourse;
    for (int id : customers) {
        assert(id >= 1 && id < instance.n);
        Node v = instance.g.nodeFromId(id);
        sortedRecourse.push_back(instance.getEdgeRecourseCost(v));
    }
    std::sort(sortedRecourse.begin(), sortedRecourse.end());

    // Create accumulated recourse vector.
    std::vector<double> accRecourse(sortedRecourse.size(), 0.0);
    accRecourse[0] = sortedRecourse[0];
    for (size_t i = 1; i < sortedRecourse.size(); i++) {
        accRecourse[i] = accRecourse[i - 1] + sortedRecourse[i];
    }

    // Compute lower bound.
    double lowerBound = 0.0;
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        // Get sum of scenario demands.
        double scenarioDemand = 0.0;
        for (int id : customers) {
            scenarioDemand += instance.scenariosMatrix[id][scenarioId];
        }

        // Get number of failures in the scenario.
        int nFailures =
            std::ceil(scenarioDemand / instance.capacity) - nVehicles;

        if (nFailures >= 1) {
            assert(static_cast<size_t>(nFailures - 1) < sortedRecourse.size());
            lowerBound += accRecourse[nFailures - 1];
        }
    }

    return lowerBound / static_cast<double>(instance.nScenarios);
}

double OptimalRecourseHelper::getRecourseCost(
    const std::vector<std::vector<int>> &routes) const {
    double totalCost = 0.0;

    for (size_t i = 0; i < routes.size(); i++) {
        PartialRoute partialRoute(instance);
        partialRoute.buildFromRoute(routes[i]);
        totalCost += getPartialRouteRecourseCost(partialRoute);
    }

    return totalCost;
}

double OptimalRecourseHelper::getPartialRouteRecourseCost(
    const PartialRoute &partialRoute) const {
    std::vector<std::vector<double>> betaDuals(
        instance.nScenarios, std::vector<double>(instance.n, 0.0));
    std::unordered_map<std::tuple<int, int, int>, double,
                       boost::hash<std::tuple<int, int, int>>>
        alphaDuals;
    return getPartialRouteRecourseCost(partialRoute, betaDuals, alphaDuals);
}

double OptimalRecourseHelper::getPartialRouteRecourseCost(
    const PartialRoute &partialRoute,
    std::vector<std::vector<double>> &betaDuals,
    std::unordered_map<std::tuple<int, int, int>, double,
                       boost::hash<std::tuple<int, int, int>>> &alphaDuals)
    const {
    if (partialRoute.entries.empty()) {
        return 0.0;
    }

    // auto started = chrono::high_resolution_clock::now();
    // double total1 = 0.0;
    // for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++)
    // {
    //     double DPrecourseCost = getPartialRouteRecourseCostInScenarioWithDP(
    //         partialRoute, scenarioId);
    //     total1 += DPrecourseCost;
    // }
    // auto done = chrono::high_resolution_clock::now();
    // double time1 =
    //     static_cast<double>(
    //         chrono::duration_cast<chrono::microseconds>(done - started)
    //             .count()) /
    //     1e6;
    // totalDPTime += time1;

    auto started = chrono::high_resolution_clock::now();
    double total2 =
        getPartialRouteRecourseCostWithLP(partialRoute, betaDuals, alphaDuals);
    auto done = chrono::high_resolution_clock::now();
    double time2 =
        static_cast<double>(
            chrono::duration_cast<chrono::microseconds>(done - started)
                .count()) /
        1e6;
    totalLPTime += time2;
    // assert(std::abs(total1 - total2) <= 1e-6);

    return total2 / static_cast<double>(instance.nScenarios);
}

double OptimalRecourseHelper::getPartialRouteRecourseCostInScenarioWithDP(
    const PartialRoute &partialRoute, int scenarioId) const {
    assert(scenarioId >= 0 && scenarioId < instance.nScenarios);
    assert(static_cast<int>(partialRoute.entries.size()) >= 1);

    // Ignore scenarios with low demands.
    double scenarioDemand = 0.0;
    for (const PartialRouteEntry &entry : partialRoute.entries) {
        for (const int id : entry.vertices) {
            assert(id >= 0 && id < instance.n);
            scenarioDemand += instance.scenariosMatrix[id][scenarioId];
        }
    }
    if (scenarioDemand <= instance.capacity) {
        return 0.0;
    }

    // Each state stores the index and the accumulated demand (after recourse).
    double bestRecourseCost = DBL_MAX;
    std::queue<std::pair<int, int>> queue;
    std::unordered_map<std::pair<int, int>, double,
                       boost::hash<std::pair<int, int>>>
        stateCost;
    queue.push({-1, 0});
    stateCost[{-1, 0}] = 0.0;

    while (queue.size() > 0) {
        std::pair<int, int> currState = queue.front();
        auto [i, demand] = currState;
        assert(stateCost.find(currState) != stateCost.end());
        assert(i >= -1 && i < static_cast<int>(partialRoute.entries.size()));
        assert(demand >= 0 && demand <= static_cast<int>(instance.capacity));

        double currRecourseCost = stateCost[currState];
        queue.pop();
        int j = i + 1;
        if (j == static_cast<int>(partialRoute.entries.size())) {
            if (currRecourseCost <= bestRecourseCost - 1e-6) {
                bestRecourseCost = currRecourseCost;
            }
            continue;
        }

        // Compute next demand.
        assert(j >= 0 && j < static_cast<int>(partialRoute.entries.size()));
        int accDemand =
            demand + partialRoute.entries[j].scenarioDemands[scenarioId];

        // Create a state that does no recourse.
        if (accDemand <= static_cast<int>(instance.capacity)) {
            std::pair<int, int> state1 = {j, accDemand};
            if (stateCost.find(state1) == stateCost.end()) {
                queue.push(state1);
                stateCost[state1] = currRecourseCost;
            } else {
                stateCost[state1] =
                    std::min(stateCost[state1], currRecourseCost);
            }
        }

        // Try to recourse at each node in the partial route set.
        // (Note: the nodes should be sorted by their recourse costs.)
        int countFailures = 0;
        int accRecourse = 0.0;
        for (const int id : partialRoute.entries[j].vertices) {
            Node v = instance.g.nodeFromId(id);
            assert(v != INVALID);
            countFailures++;
            accRecourse += instance.getEdgeRecourseCost(v);
            std::pair<int, int> state2 = {
                j, std::max(accDemand - countFailures *
                                            static_cast<int>(instance.capacity),
                            0)};

            if (state2.second > static_cast<int>(instance.capacity)) {
                continue;
            }

            if (currRecourseCost + accRecourse <= bestRecourseCost) {
                if (stateCost.find(state2) == stateCost.end()) {
                    queue.push(state2);
                    stateCost[state2] = currRecourseCost + accRecourse;
                } else {
                    stateCost[state2] = std::min(
                        stateCost[state2], currRecourseCost + accRecourse);
                }
            }

            // No need to examine more nodes.
            if (state2.second == 0) {
                break;
            }
        }
    }

    stateCost.clear();
    return bestRecourseCost;
}

double OptimalRecourseHelper::getPartialRouteRecourseCostWithLP(
    const PartialRoute &partialRoute,
    std::vector<std::vector<double>> &betaDuals,
    std::unordered_map<std::tuple<int, int, int>, double,
                       boost::hash<std::tuple<int, int, int>>> &alphaDuals)
    const {
    GRBModel model(env);
    model.set(GRB_StringAttr_ModelName, "OptimalRecourseModel");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model.set(GRB_DoubleParam_TimeLimit, 300);
    model.set(GRB_IntParam_Threads, 1);
    // model.set(GRB_IntParam_Aggregate, 2);

    // Precompute sum of route demands.
    bool noFailures = true;
    // std::vector<std::vector<int>> sumDemands(
    //     partialRoute.entries.size(), std::vector<int>(instance.nScenarios,
    //     0));
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        if (partialRoute.totalScenarioDemands[scenarioId] <=
            instance.capacity) {
            continue;
        }

        noFailures = false;
        // for (size_t i = 0; i < partialRoute.entries.size(); i++) {
        //     int accLoad = (i == 0) ? 0 : sumDemands[i - 1][scenarioId];
        //     accLoad += partialRoute.entries[i].scenarioDemands[scenarioId];
        //     sumDemands[i][scenarioId] = accLoad;
        // }
    }

    if (noFailures) {
        return 0.0;
    }

    // Customer variables.
    std::vector<std::vector<GRBVar>> customerVar(
        instance.n, std::vector<GRBVar>(instance.nScenarios));
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        if (partialRoute.totalScenarioDemands[scenarioId] <=
            instance.capacity) {
            continue;
        }

        for (size_t i = 0; i < partialRoute.entries.size(); i++) {
            for (int id : partialRoute.entries[i].vertices) {
                assert(id >= 1 && id < instance.n);
                Node v = instance.g.nodeFromId(id);
                customerVar[id][scenarioId] = model.addVar(
                    0.0, 1.0, instance.getEdgeRecourseCost(v), GRB_CONTINUOUS,
                    "customer_recourse_" + std::to_string(id) + "_scen_" +
                        std::to_string(scenarioId));
            }
        }
    }
    model.update();

    // Add constraints.
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        if (partialRoute.totalScenarioDemands[scenarioId] <=
            instance.capacity) {
            continue;
        }

        int limit = std::ceil(partialRoute.totalScenarioDemands[scenarioId] /
                              instance.capacity);

        for (int t = 2; t <= limit; t++) {
            // Sliding window that tracks subpartial routes.
            size_t i = 0;
            size_t j = 0;
            double currDemand = 0.0;

            while (i < partialRoute.entries.size()) {
                // Find the end index of the subpartial route.
                while (j < partialRoute.entries.size() &&
                       std::ceil(currDemand / instance.capacity) < t) {
                    currDemand +=
                        partialRoute.entries[j].scenarioDemands[scenarioId];
                    j++;
                }

                // Increase the start index.
                while (
                    std::ceil(
                        (currDemand -
                         partialRoute.entries[i].scenarioDemands[scenarioId]) /
                        instance.capacity) >= t) {
                    currDemand -=
                        partialRoute.entries[i].scenarioDemands[scenarioId];
                    i++;
                }

                assert(i < j);

                // Found a minimal set, add constraint only if `ceil(currDemand
                // / C) = t`.
                if (static_cast<int>(
                        std::ceil(currDemand / instance.capacity)) == t) {
                    GRBLinExpr expr = 0.0;
                    for (size_t k = i; k < j; k++) {
                        for (int id : partialRoute.entries[k].vertices) {
                            expr += customerVar[id][scenarioId];
                        }
                    }
                    std::string consName =
                        "cons_" + std::to_string(scenarioId) + "," +
                        std::to_string(i) + "," + std::to_string(j);
                    model.addConstr(expr >= t - 1, consName);
                }

                currDemand -=
                    partialRoute.entries[i].scenarioDemands[scenarioId];
                i++;
            }
        }

        // Different way of adding constraints (slower).
        // for (size_t i = 0; i < partialRoute.entries.size(); i++) {
        //     // We use this vector to avoid adding redundant constraints.
        //     std::vector<bool> addedConstraint(failLimit + 1, false);
        //     addedConstraint[0] = true;
        //     addedConstraint[1] = true;

        //     for (size_t j = i; j < partialRoute.entries.size(); j++) {
        //         int prev = (i == 0) ? 0 : sumDemands[i - 1][scenarioId];
        //         int subrouteDemand = sumDemands[j][scenarioId] - prev;
        //         int currentFail = std::ceil(subrouteDemand /
        //         instance.capacity);

        //         if (addedConstraint[currentFail]) {
        //             continue;
        //         }
        //         addedConstraint[currentFail] = true;

        //         GRBLinExpr expr = 0.0;
        //         for (size_t k = i; k <= j; k++) {
        //             for (int id : partialRoute.entries[k].vertices) {
        //                 expr += customerVar[id][scenarioId];
        //             }
        //         }
        //         std::string consName = "cons_" + std::to_string(scenarioId) +
        //                                "," + std::to_string(i) + "," +
        //                                std::to_string(j);
        //         model.addConstr(
        //             expr >= std::ceil(subrouteDemand / instance.capacity) -
        //             1, consName);
        //     }
        // }
    }

    model.update();
    model.optimize();

    if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
        std::cout << "Error when computing optimal recourse cost." << std::endl;
        exit(1);
    }

    double primalObj = model.get(GRB_DoubleAttr_ObjVal);

    // Now we need to construct a dual solution with minimum support.
    // First, we get the beta coefs.
    double dualObj = 0.0;
    std::vector<std::vector<bool>> seletedBetas(std::vector<std::vector<bool>>(
        instance.nScenarios,
        std::vector<bool>(instance.n, false))); // This will help later.
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        if (partialRoute.totalScenarioDemands[scenarioId] <=
            instance.capacity) {
            continue;
        }

        for (size_t i = 0; i < partialRoute.entries.size(); i++) {
            for (int id : partialRoute.entries[i].vertices) {
                // For bound constraints, the dual is in the reduced cost.
                betaDuals[scenarioId][id] = std::min(
                    customerVar[id][scenarioId].get(GRB_DoubleAttr_RC), 0.0);
                dualObj += betaDuals[scenarioId][id];

                if (betaDuals[scenarioId][id] <= -1e-5) {
                    seletedBetas[scenarioId][id] = true;
                }
            }
        }
    }

    // Next, we get the alpha coefs.
    GRBConstr *modelCons = 0;
    modelCons = model.getConstrs();
    double nuCoef = 0.0;
    for (int consIdx = 0; consIdx < model.get(GRB_IntAttr_NumConstrs);
         ++consIdx) {
        double alpha = modelCons[consIdx].get(GRB_DoubleAttr_Pi);
        double rhs = modelCons[consIdx].get(GRB_DoubleAttr_RHS);

        if (alpha >= 1e-5) {
            std::string consName =
                modelCons[consIdx].get(GRB_StringAttr_ConstrName);
            int pos = consName.find_first_of('_');
            std::stringstream ss(consName.substr(pos + 1));

            int scenarioId, startIdx, endIdx;
            char comma = ',';
            ss >> scenarioId >> comma >> startIdx >> comma >> endIdx;

            // Check if we need to update the duals (this is tricky, see the
            // paper).
            double supportSize = 0.0;
            for (size_t i = startIdx; i < endIdx; i++) {
                for (int id : partialRoute.entries[i].vertices) {
                    if (seletedBetas[scenarioId][id]) {
                        supportSize += 1.0;
                    }
                }
            }

            if (rhs - supportSize <= 0.0) {
                assert(std::abs(rhs - supportSize) <= 1e-6);

                // Compute step size.
                double stepSize = alpha;
                for (size_t i = startIdx; i < endIdx; i++) {
                    for (int id : partialRoute.entries[i].vertices) {
                        if (betaDuals[scenarioId][id] <= -1e-5) {
                            stepSize =
                                std::min(stepSize, -betaDuals[scenarioId][id]);
                        }
                    }
                }

                // Update beta coefs and support.
                alpha -= stepSize;
                for (size_t i = startIdx; i < endIdx; i++) {
                    for (int id : partialRoute.entries[i].vertices) {
                        if (betaDuals[scenarioId][id] <= -1e-5) {
                            betaDuals[scenarioId][id] += stepSize;

                            if (std::abs(betaDuals[scenarioId][id]) <= 1e-5) {
                                seletedBetas[scenarioId][id] = false;
                            }
                        }
                    }
                }
            }

            auto t = std::make_tuple(scenarioId, startIdx, endIdx);
            if (alpha >= 1e-6) {
                alphaDuals[t] = alpha;
            }

            // Update dual obj.
            dualObj += alpha * rhs;
        }
    }
    delete[] modelCons;

    assert(std::abs(dualObj - primalObj) <= 1e-6);
    return primalObj;
}
