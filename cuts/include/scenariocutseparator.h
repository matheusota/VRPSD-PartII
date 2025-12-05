#ifndef SCENARIOCUTSEPARATOR_H
#define SCENARIOCUTSEPARATOR_H

// #include <gurobi_c++.h>
#include "mygraphlib.h"
#include "svrpinstance.h"
#include "params.h"
#include <cassert>
#include <list>
#include <utility>
#include <ctime> // For CLOCKS_PER_SEC
#include <chrono>
#include <unordered_set>
#include <algorithm>
#include "cutdata.h"

class ScenarioCutSeparator {
  public:
    int nCuts = 0;
    double time = 0.0;

    ScenarioCutSeparator(const SVRPInstance &instance, const Params &params);

    // Return the number of separated cuts.
    int separateScenarioCut(int scenarioId,
                            const std::vector<std::vector<int>> &customerSets,
                            const EdgeValueMap &xValue, double scenarioValue,
                            std::vector<CutData> &separatedCuts);

  private:
    const SVRPInstance &instance;
    const Params &params;
    double EpsForIntegrality = 1e-4;

    void addCustomerSetToCut(const EdgeValueMap &xValue,
                             const std::vector<int> customers, int nVehicles,
                             double coef, CutData &cutData);

    int addCut(const EdgeValueMap &xValue, const NodeValueMap &yValue,
               int nScenariosCut, double RHS,
               const std::unordered_set<int> &customerSet,
               std::vector<CutData> &separatedCuts);
};

#endif // SCENARIOCUTSEPARATOR_H