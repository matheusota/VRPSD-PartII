#include "scenariocutseparator.h"

ScenarioCutSeparator::ScenarioCutSeparator(const SVRPInstance &instance,
                                           const Params &params)
    : instance(instance), params(params) {}

// Here we add the term coef * (k(S) - x(\delta(S)) / 2) to the LHS of the
// inequality. In fact, since 2 * x(E(S)) + x(\delta(S)) = 2 * |S|,
// we add the term coef * (k(S) - |S| + x(E(S))) instead.
void ScenarioCutSeparator::addCustomerSetToCut(const EdgeValueMap &xValue,
                                               const std::vector<int> customers,
                                               int nVehicles, double coef,
                                               CutData &cutData) {
    cutData.RHS += coef * (nVehicles - static_cast<double>(customers.size()));

    for (size_t i = 0; i < customers.size(); i++) {
        for (size_t j = i + 1; j < customers.size(); j++) {
            Node u = instance.g.nodeFromId(customers[i]);
            Node v = instance.g.nodeFromId(customers[j]);
            Edge e = findEdge(instance.g, u, v);
            cutData.edgePairs.push_back(std::make_pair(e, -coef));
            cutData.LHS += -coef * xValue[e];
        }
    }
}

int ScenarioCutSeparator::separateScenarioCut(
    int scenarioId, const std::vector<std::vector<int>> &customerSets,
    const EdgeValueMap &xValue, double scenarioValue,
    std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    // First we order the customer sets by their lowest recourse values.
    std::vector<std::tuple<int, double, int, int>> orderedCustomerSets;
    for (size_t setId = 0; setId < customerSets.size(); setId++) {
        const std::vector<int> &customers = customerSets[setId];

        // Get total demand and lowest recourse.
        double lowestRecourse = 1e6;
        double totalScenarioDemand = 0;
        double totalDemand = 0;

        for (const int id : customers) {
            Node v = instance.g.nodeFromId(id);
            lowestRecourse =
                std::min(lowestRecourse, instance.getEdgeRecourseCost(v));
            totalScenarioDemand += instance.scenariosMatrix[id][scenarioId];
            totalDemand += instance.demand[v];
        }

        int scenarioNVehicles = static_cast<int>(
            std::ceil(totalScenarioDemand / instance.capacity));
        int nVehicles =
            static_cast<int>(std::ceil(totalDemand / instance.capacity));

        if (scenarioNVehicles <= nVehicles) {
            continue;
        }

        orderedCustomerSets.push_back(std::make_tuple(
            customers.size(), -lowestRecourse, scenarioNVehicles, setId));
    }
    std::sort(orderedCustomerSets.begin(), orderedCustomerSets.end());

    if (orderedCustomerSets.empty()) {
        return 0;
    }

    // Now we construct the corresponding scenario cut.
    separatedCuts.emplace_back(CutData());
    CutData &cutData = separatedCuts.back();
    cutData.RHS = 0.0;
    cutData.LHS = scenarioValue;
    cutData.sense = '>';
    cutData.name = "ScenarioCut " + std::to_string(scenarioId);
    cutData.scenarioPairs.push_back(std::make_pair(scenarioId, 1.0));

    std::vector<double> customerLHS(instance.n, 0.0);
    for (const std::tuple<int, double, int, int> &t : orderedCustomerSets) {
        double coef = 1e6;
        int newNVehicles = std::get<2>(t);
        int setId = std::get<3>(t);

        for (const int id : customerSets[setId]) {
            coef = std::min(
                coef, instance.getEdgeRecourseCost(instance.g.nodeFromId(id)) -
                          customerLHS[id]);
        }

        if (coef <= EpsForIntegrality) {
            continue;
        }

        addCustomerSetToCut(xValue, customerSets[setId], newNVehicles, coef,
                            cutData);

        for (const int id : customerSets[setId]) {
            customerLHS[id] += coef;
        }
    }

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    if (cutData.LHS >= cutData.RHS - EpsForIntegrality) {
        separatedCuts.pop_back();
        return 0;
    } else {
        cutData.violation = cutData.RHS - cutData.LHS;
        nCuts++;
        return 1;
    }
}