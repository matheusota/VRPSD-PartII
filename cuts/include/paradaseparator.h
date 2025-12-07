#ifndef PARADASEPARATOR_H
#define PARADASEPARATOR_H

// #include <gurobi_c++.h>
#include "mygraphlib.h"
#include "svrpinstance.h"
#include "params.h"
#include <cassert>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <list>
#include <utility>
#include <ctime> // For CLOCKS_PER_SEC
#include <chrono>
#include "cutdata.h"

class ParadaSeparator {
  public:
    int nSetCuts = 0;
    int nPathCuts = 0;
    double time = 0.0;

    ParadaSeparator(const SVRPInstance &instance, const Params &params);

    // These functions return the number of separated cuts.
    int separateSetCut(const EdgeValueMap &xValue,
                       const NodeValueMap &recourseValue,
                       const std::vector<int> &customers,
                       std::vector<CutData> &separatedCuts);

  private:
    const SVRPInstance &instance;
    const Params &params;
    double EpsForIntegrality = 1e-4;

    int addSetCut(const EdgeValueMap &xValue, const NodeValueMap &recourseValue,
                  const std::vector<int> &customers,
                  std::vector<CutData> &separatedCuts);
};

#endif // PARADASEPARATOR_H