#ifndef PARTIALROUTESEPARATOR_H
#define PARTIALROUTESEPARATOR_H

#include "mygraphlib.h"
#include "svrpinstance.h"
#include "graphviewer.h"
#include <cassert>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/dfs.h>
#include <list>
#include <utility>
#include <set>
#include <algorithm>
#include <random>
#include <chrono>
#include "params.h"
#include "partialroute.h"
#include "cutdata.h"

typedef Dfs<Graph> DFS;

class PartialRouteSeparator {
  public:
    int nCuts = 0;
    int elementaryCuts = 0;
    double time = 0.0;

    PartialRouteSeparator(const SVRPInstance &instance, const Params &params);

    int separatePartialRouteCuts(const EdgeValueMap &xValue,
                                 const NodeValueMap &recourseValue,
                                 const std::vector<PartialRoute> &partialRoutes,
                                 std::vector<CutData> &separatedCuts);

  private:
    const SVRPInstance &instance;
    const Params &params;
    double EpsForIntegrality = 1e-4;

    void setCutCoefficientsNew(const PartialRoute &partialRoute,
                               double recourseCost, EdgeValueMap &edgeCoefs,
                               double &RHS, bool simpleAdherence,
                               bool improved);
    // Set coefficients x(E(S)).
    void setCoefficientsInsideSet(const std::vector<int> &vertices, double coef,
                                  EdgeValueMap &edgeCoefs);
    // Set coefficients x(v, S).
    void setCoefficientsToSet(int from, const std::vector<int> &vertices,
                              double coef, EdgeValueMap &edgeCoefs);
    bool addCutFromCoefs(const EdgeValueMap &xValue,
                         const std::vector<Node> &nodesToConsider,
                         const NodeValueMap &recourseValue,
                         const EdgeValueMap &edgeCoefs, double RHS,
                         std::vector<CutData> &separatedCuts);
};

#endif // GRBPARTIALROUTESEPARATOR_H
