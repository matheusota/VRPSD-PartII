#include "optimalrecoursehelper.h"
#include "svrpinstance.h"

OptimalRecourseHelper::OptimalRecourseHelper(const SVRPInstance &instance_)
    : instance(instance_) {}

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
        totalCost += getRouteRecourseCost(routes[i]);
    }

    return totalCost;
}

double OptimalRecourseHelper::getRouteRecourseCost(
    const std::vector<int> &route) const {
    double total = 0.0;

    if (route.empty()) {
        return 0.0;
    }

    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        double DPrecourseCost =
            getRouteRecourseCostInScenarioWithDP(route, scenarioId);
        // std::vector<std::unordered_set<int>> tmp;
        // double LPrecourseCost =
        //     getRouteRecourseCostInScenarioWithLP(route, scenarioId, tmp);
        // assert(std::abs(LPrecourseCost - DPrecourseCost) <= 1e-6);
        total += DPrecourseCost;
    }

    return total / static_cast<double>(instance.nScenarios);
}

double OptimalRecourseHelper::getRouteRecourseCostInScenarioWithDP(
    const std::vector<int> &route, int scenarioId) const {
    // Ignore scenarios with low demands.
    double scenarioDemand = 0.0;
    for (const int id : route) {
        assert(id >= 0 && id < instance.n);
        scenarioDemand += instance.scenariosMatrix[id][scenarioId];
    }
    if (scenarioDemand <= instance.capacity) {
        return 0.0;
    }

    assert(scenarioId >= 0 && scenarioId < instance.nScenarios);
    assert(static_cast<int>(route.size()) >= 1);

    double bestRecourseCost =
        instance.classicalRecourseHelper.getRouteRecourseCostInScenario(
            route, scenarioId, false);
    std::queue<std::pair<int, int>> queue;
    std::unordered_map<std::pair<int, int>, double,
                       boost::hash<std::pair<int, int>>>
        stateCost;
    std::unordered_map<std::pair<int, int>, int,
                       boost::hash<std::pair<int, int>>>
        stateFailures;
    queue.push({-1, 0});
    stateCost[{-1, 0}] = 0.0;
    stateFailures[{-1, 0}] = 0;

    while (queue.size() > 0) {
        std::pair<int, int> currState = queue.front();
        auto [i, demand] = currState;
        assert(stateCost.find(currState) != stateCost.end() &&
               stateFailures.find(currState) != stateFailures.end());
        assert(i >= -1 && i < static_cast<int>(route.size()));
        assert(demand >= 0 && demand <= static_cast<int>(instance.capacity));

        double currRecourseCost = stateCost[currState];
        int currFailures = stateFailures[currState];
        queue.pop();
        int j = i + 1;
        if (j == static_cast<int>(route.size())) {
            if (currRecourseCost <= bestRecourseCost - 1e-6) {
                bestRecourseCost = currRecourseCost;
            }
            continue;
        }

        // For each state, we have two decisions, return to the depot or
        // not.
        assert(route[j] >= 0 && route[j] < instance.n);
        int nextDemand =
            demand + instance.scenariosMatrix[route[j]][scenarioId];
        std::pair<int, int> state1 = {j, nextDemand};
        std::pair<int, int> state2 = {
            j, std::max(nextDemand - static_cast<int>(instance.capacity), 0)};
        double nextRecourseCost =
            currRecourseCost +
            instance.getEdgeRecourseCost(instance.g.nodeFromId(route[j]));
        int nextFailures = currFailures + 1;

        if (nextDemand <= static_cast<int>(instance.capacity)) {
            if (stateCost.find(state1) == stateCost.end()) {
                queue.push(state1);
                stateCost[state1] = currRecourseCost;
                stateFailures[state1] = currFailures;
            } else {
                stateCost[state1] =
                    std::min(stateCost[state1], currRecourseCost);
                stateFailures[state1] =
                    std::min(stateFailures[state1], currFailures);
            }
        }

        if (nextDemand >= 1 && nextRecourseCost <= bestRecourseCost) {
            if (stateCost.find(state2) == stateCost.end()) {
                queue.push(state2);
                stateCost[state2] = nextRecourseCost;
                stateFailures[state2] = nextFailures;
            } else {
                stateCost[state2] =
                    std::min(stateCost[state2], nextRecourseCost);
                stateFailures[state2] =
                    std::min(stateFailures[state2], nextFailures);
            }
        }
    }

    stateCost.clear();
    return bestRecourseCost;
}

double OptimalRecourseHelper::getPartialRouteRecourseCost(
    const PartialRoute &partialRoute) const {
    if (partialRoute.entries.empty()) {
        return 0.0;
    }

    double total = 0.0;
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        double DPrecourseCost = getPartialRouteRecourseCostInScenarioWithDP(
            partialRoute, scenarioId);
        total += DPrecourseCost;
    }

    return total / static_cast<double>(instance.nScenarios);
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
        int nextDemand = 0;
        assert(j >= 0 && j < static_cast<int>(partialRoute.entries.size()));
        for (const int id : partialRoute.entries[j].vertices) {
            nextDemand += instance.scenariosMatrix[id][scenarioId];
        }
        int accDemand = demand + nextDemand;

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

double OptimalRecourseHelper::getRouteRecourseCostInScenarioWithLP(
    const std::vector<int> &route, int scenarioId,
    std::vector<std::unordered_set<int>> &tightSubroutes) const {
    GRBEnv env(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    GRBModel model(env);

    model.set(GRB_StringAttr_ModelName, "OptimalRecourseModel");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model.set(GRB_DoubleParam_TimeLimit, 300);
    model.set(GRB_DoubleParam_MIPGap, 1e-9);
    model.set(GRB_IntParam_Threads, 1);
    env.set(GRB_IntParam_LogToConsole, 0);

    // Customer variables.
    std::vector<GRBVar> customerVar(route.size());
    for (size_t i = 0; i < route.size(); i++) {
        Node v = instance.g.nodeFromId(route[i]);

        customerVar[i] = model.addVar(
            0.0, GRB_INFINITY, instance.getEdgeRecourseCost(v), GRB_CONTINUOUS,
            "customer_recourse_" + std::to_string(route[i]) + "_scen_" +
                std::to_string(scenarioId));
    }
    model.update();

    // Precompute sum of route demands.
    std::vector<int> sumDemands(route.size(), 0);
    for (size_t i = 0; i < route.size(); i++) {
        int prev = (i == 0) ? 0 : sumDemands[i - 1];
        assert(route[i] >= 0 && route[i] < instance.n);
        sumDemands[i] = prev + instance.scenariosMatrix[route[i]][scenarioId];
    }

    // Subroute constraints.
    for (size_t i = 0; i < route.size(); i++) {
        for (size_t j = i; j < route.size(); j++) {
            int prev = (i == 0) ? 0 : sumDemands[i - 1];
            int subrouteDemand = sumDemands[j] - prev;

            if (subrouteDemand <= static_cast<int>(instance.capacity)) {
                continue;
            }

            GRBLinExpr expr = 0.0;
            std::string consName = "route_";
            for (size_t k = i; k <= j; k++) {
                expr += customerVar[k];
                consName += std::to_string(route[k]) + ",";
            }
            consName += "#" + std::to_string(i) + "," + std::to_string(j);

            model.addConstr(
                expr >= std::ceil(subrouteDemand / instance.capacity) - 1,
                consName);
        }
    }

    model.update();
    model.optimize();

    if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
        std::cout << "Error when computing optimal recourse cost." << std::endl;
        exit(1);
    }

    double recourseCost = model.get(GRB_DoubleAttr_ObjVal);

    GRBConstr *modelCons = 0;
    modelCons = model.getConstrs();

    for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); ++i) {
        if (std::abs(modelCons[i].get(GRB_DoubleAttr_Pi)) > 1e-6) {
            std::string consName = modelCons[i].get(GRB_StringAttr_ConstrName);
            int pos = consName.find_first_of('#');
            std::string tmp = consName.substr(pos + 1);
            pos = tmp.find_first_of(',');
            int first = std::stoi(tmp.substr(0, pos));
            int second = std::stoi(tmp.substr(pos + 1));

            std::unordered_set<int> subroute;
            for (int j = first; j <= second; j++) {
                subroute.insert(route[j]);
            }
            tightSubroutes.push_back(subroute);
        }
    }

    return recourseCost;
}