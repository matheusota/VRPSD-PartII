#ifndef BASICCALLBACK_H
#define BASICCALLBACK_H

#include <gurobi_c++.h>
#include "mygraphlib.h"
#include "svrpinstance.h"
#include "params.h"
#include "graphviewer.h"
#include <cassert>
#include <lemon/list_graph.h>
#include <list>
#include <utility>
#include <ctime> // For CLOCKS_PER_SEC
#include <chrono>
#include "cvrpsepseparator.h"
#include "gendreauseparator.h"
#include "cutdata.h"

class BasicCallback : public GRBCallback {
  public:
    double rootBound;
    std::vector<std::tuple<GRBLinExpr, char, double, std::string>> lazyVector;

    BasicCallback(SVRPInstance &instance, Params &params, GRBModel &model,
                  EdgeGRBVarMap &x, std::vector<GRBVar> &scenariosVars,
                  CVRPSEPSeparator *cvrpsepSeparator,
                  GendreauSeparator *gendreauSeparator);
    void callback();

  private:
    double (GRBCallback::*solutionValue)(GRBVar);

    SVRPInstance &instance;
    Params &params;
    GRBModel &model;
    EdgeGRBVarMap &x;
    std::vector<GRBVar> &scenariosVars;
    CVRPSEPSeparator *cvrpsepSeparator;
    GendreauSeparator *gendreauSeparator;
    const int noImprovLimit = 5;
    const double noImprovEps = 1e-4;
    double prevDualBound;
    int prevNodeCount;
    int countNoImprov;

    void updateImprovData(int nodeCount, double dualBound);
};

#endif // BASICCALLBACK_H
