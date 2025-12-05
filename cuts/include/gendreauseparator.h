#ifndef GENDREAUSEPARATOR_H
#define GENDREAUSEPARATOR_H

// #include <gurobi_c++.h>
#include "mygraphlib.h"
#include "svrpinstance.h"
#include "params.h"
#include <cassert>
#include <lemon/list_graph.h>
#include <lemon/maps.h>
#include <list>
#include <utility>
#include <ctime> // For CLOCKS_PER_SEC
#include <chrono>
#include "cutdata.h"

class GendreauSeparator {
  public:
    int nCuts = 0;
    double time = 0.0;

    GendreauSeparator(const SVRPInstance &instance, const Params &params);

    int separateGendreauCutsPerScenario(
        const EdgeValueMap &xValue,
        const std::vector<double> &scenariosRecourse,
        std::vector<CutData> &separatedCuts);

    int separateAggregatedGendreauCut(const EdgeValueMap &xValue,
                                      double totalRecourse,
                                      std::vector<CutData> &separatedCuts);

  private:
    const SVRPInstance &instance;
    const Params &params;
    double EpsForIntegrality = 1e-4;

    bool captureRoutes(const EdgeValueMap &xValue,
                       std::vector<std::vector<int>> &routes);
    int addGendreauCutForScenario(const EdgeValueMap &xValue,
                                  const std::vector<double> &scenariosRecourse,
                                  const std::vector<std::vector<int>> &routes,
                                  int scenarioId,
                                  std::vector<CutData> &separatedCuts);

    int addAggregatedGendreauCut(const EdgeValueMap &xValue,
                                 double totalRecourse,
                                 const std::vector<std::vector<int>> &routes,
                                 std::vector<CutData> &separatedCuts);
};

#endif // GENDREAUSEPARATOR_H