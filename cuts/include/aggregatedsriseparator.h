#ifndef AGGREGATEDSRISEPARATOR_H
#define AGGREGATEDSRISEPARATOR_H

#include <gurobi_c++.h>
#include "mygraphlib.h"
#include "svrpinstance.h"
#include "params.h"
#include <cassert>
#include <list>
#include <utility>
#include <ctime> // For CLOCKS_PER_SEC
#include <chrono>
#include <unordered_set>
#include <lemon/connectivity.h>
#include "cutdata.h"
#include "svrpsolution.h"
#include "partialroutesbuilder.h"

class AggregatedSRISeparator {
  public:
    int nCuts = 0;
    double time = 0.0;

    AggregatedSRISeparator(const SVRPInstance &instance, const Params &params);

    int addCutFromSet(const EdgeValueMap &xValue, const NodeValueMap &yValue,
                      const std::vector<int> &customers,
                      std::vector<CutData> &separatedCuts);
    int flowHeuristicSeparation(const EdgeValueMap &xValue,
                                const NodeValueMap &yValue,
                                std::vector<CutData> &separatedCuts);
    void addCriticalScenarios(const std::vector<int> scenarioIds) {
        criticalScenarios.insert(scenarioIds);
    }

  private:
    const SVRPInstance &instance;
    const Params &params;
    double EpsForIntegrality = 1e-4;
    OptimalRecourseHelper optimalRecourseHelper;

    // Digraph for the network flow heuristic.
    double flowRoundFactor = 1e3;
    Digraph flowGraph;
    DNode newSource;

    // Use this to store critical scenarios when calling flow separation.
    std::unordered_set<std::vector<int>, hash_vector> criticalScenarios;

    int flowSeparationFromScenarios(const EdgeValueMap &xValue,
                                    const NodeValueMap &yValue,
                                    const std::vector<int> &scenarios,
                                    std::vector<CutData> &separatedCuts);
    int addCutFromSetSimple(const EdgeValueMap &xValue,
                            const NodeValueMap &yValue,
                            const std::vector<int> &customers,
                            std::vector<CutData> &separatedCuts);
    int addCutFromSetWithDual(const EdgeValueMap &xValue,
                              const NodeValueMap &yValue,
                              const std::vector<int> &customers,
                              std::vector<CutData> &separatedCuts);
};

#endif // AGGREGATEDSRISEPARATOR_H