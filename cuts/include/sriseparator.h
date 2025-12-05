#ifndef SRISEPARATOR_H
#define SRISEPARATOR_H

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
#include "cutdata.h"
#include "cvrpsepseparator.h"
#include "svrpsolution.h"

enum SRI_SEPARATION_STATUS { FOUND_CUTS, NO_CUTS, FAILED };

class SRISeparator {
  public:
    int nCuts = 0;
    double time = 0.0;

    SRISeparator(const SVRPInstance &instance_, const Params &params_,
                 CVRPSEPSeparator &cvrpsepSeparator);

    // Return the number of separated cuts.
    int flowSeparation(const EdgeValueMap &xValue,
                       const NodeVectorValueMap &yValue,
                       std::vector<CutData> &separatedCuts);
    int heuristicSeparation(const EdgeValueMap &xValue,
                            const NodeVectorValueMap &yValue,
                            std::vector<CutData> &separatedCuts);
    int mipSeparation(const EdgeValueMap &xValue,
                      const NodeVectorValueMap &yValue,
                      std::vector<CutData> &separatedCuts,
                      SRI_SEPARATION_STATUS &status);
    int integerSeparation(const EdgeValueMap &xValue,
                          const NodeVectorValueMap &yValue,
                          std::vector<CutData> &separatedCuts);
    int separateSRIFromCustomerSet(const EdgeValueMap &xValue,
                                   const NodeVectorValueMap &yValue,
                                   const std::vector<int> &customers,
                                   std::vector<CutData> &separatedCuts);

  private:
    const SVRPInstance &instance;
    const Params &params;
    double EpsForIntegrality = 1e-4;
    CVRPSEPSeparator &cvrpsepSeparator;

    int addSRICutsFromRoute(const EdgeValueMap &xValue,
                            const NodeVectorValueMap &yValue,
                            const std::vector<int> &route,
                            std::vector<CutData> &separatedCuts);

    int addCut(const EdgeValueMap &xValue, const NodeVectorValueMap &yValue,
               const std::vector<std::pair<int, int>> &scenarioBounds,
               const std::vector<int> &customers,
               std::vector<CutData> &separatedCuts);

    int addSRIFromCustomerSet(const EdgeValueMap &xValue,
                              const NodeVectorValueMap &yValue,
                              const std::vector<int> &customers,
                              std::vector<CutData> &separatedCuts);

    SRI_SEPARATION_STATUS
    singleScenarioMIPModel(const EdgeValueMap &xValue,
                           const NodeVectorValueMap &yValue, int scenarioId,
                           std::vector<int> &customers);

    // Not used.
    int integerSeparation2(const EdgeValueMap &xValue,
                           const NodeValueMap &yValue,
                           std::vector<CutData> &separatedCuts);

    int addBendersCut(const EdgeValueMap &xValue, const NodeValueMap &yValue,
                      const std::vector<int> &route,
                      std::vector<CutData> &separatedCuts);
};

#endif // SRISEPARATOR_H