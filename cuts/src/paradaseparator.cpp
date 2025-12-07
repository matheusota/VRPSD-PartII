#include "paradaseparator.h"

ParadaSeparator::ParadaSeparator(const SVRPInstance &instance,
                                 const Params &params)
    : instance(instance), params(params) {}

int ParadaSeparator::separateSetCut(const EdgeValueMap &xValue,
                                    const NodeValueMap &recourseValue,
                                    const std::vector<int> &customers,
                                    std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    int sepCuts = addSetCut(xValue, recourseValue, customers, separatedCuts);

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return sepCuts;
}

// The return value of this indicates if the inequality is violated or not.
int ParadaSeparator::addSetCut(const EdgeValueMap &xValue,
                               const NodeValueMap &recourseValue,
                               const std::vector<int> &customers,
                               std::vector<CutData> &separatedCuts) {
    double totalDemand = 0.0;
    for (int id : customers) {
        Node v = instance.g.nodeFromId(id);
        assert(id >= 1 && id < instance.n && v != INVALID);
        totalDemand += instance.demand[v];
    }
    int nVehicles =
        (params.basicVRPSD)
            ? 1
            : static_cast<int>(std::ceil(totalDemand / instance.capacity));

    // Compute recourse.
    double lowerBound = instance.optimalRecourseHelper.getRecourseLowerBound(
        customers, nVehicles);
    if (lowerBound <= EpsForIntegrality) {
        return 0;
    }

    // Construct cut.
    separatedCuts.emplace_back(CutData());
    CutData &cutData = separatedCuts.back();
    cutData.name = "SetCut";
    cutData.RHS =
        lowerBound * (nVehicles + 1 - static_cast<int>(customers.size()));
    cutData.sense = '>';
    cutData.LHS = 0.0;

    for (int id : customers) {
        Node v = instance.g.nodeFromId(id);
        assert(v != INVALID);
        cutData.nodePairs.push_back({instance.g.id(v), 1.0});
        cutData.LHS += recourseValue[v];
        cutData.customers.push_back(id);
    }

    for (size_t j = 0; j < customers.size(); j++) {
        for (size_t k = j + 1; k < customers.size(); k++) {
            Edge e = findEdge(instance.g, instance.g.nodeFromId(customers[j]),
                              instance.g.nodeFromId(customers[k]));
            assert(e != INVALID);
            cutData.edgePairs.push_back({e, -lowerBound});
            cutData.LHS -= lowerBound * xValue[e];
        }
    }

    if (cutData.LHS >= cutData.RHS - EpsForIntegrality) {
        separatedCuts.pop_back();
        return 0;
    } else {
        cutData.violation = cutData.RHS - cutData.LHS;
        nSetCuts++;
        return 1;
    }
}