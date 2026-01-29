#ifndef PARTIALROUTESEPARATOR_H
#define PARTIALROUTESEPARATOR_H

#include "mygraphlib.h"
#include "svrpinstance.h"
#include "graphviewer.h"
#include <cassert>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <boost/functional/hash.hpp>
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

    PartialRouteSeparator(SVRPInstance &instance, const Params &params);

    int separatePartialRouteCuts(const EdgeValueMap &xValue,
                                 const NodeValueMap &recourseValue,
                                 const std::vector<PartialRoute> &partialRoutes,
                                 std::vector<CutData> &separatedCuts);

  private:
    SVRPInstance &instance;
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
                         const NodeValueMap &recourseValue,
                         const EdgeValueMap &edgeCoefs,
                         const NodeValueMap &nodeCoefs, double RHS,
                         std::vector<CutData> &separatedCuts);

    void setCutCoefficientsFromDual(
        const PartialRoute &partialRoute, const EdgeValueMap &xValue,
        double recourseCost,
        const std::unordered_map<std::pair<int, int>, double,
                                 boost::hash<std::pair<int, int>>> &betaDuals,
        const std::unordered_map<std::tuple<int, int, int>, double,
                                 boost::hash<std::tuple<int, int, int>>>
            &alphaDuals,
        EdgeValueMap &edgeCoefs, NodeValueMap &nodeCoefs, double &RHS) const;
};

#endif // GRBPARTIALROUTESEPARATOR_H
