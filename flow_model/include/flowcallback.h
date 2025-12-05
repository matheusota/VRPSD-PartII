#ifndef FLOWCALLBACK_H
#define FLOWCALLBACK_H

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
#include "cutdata.h"

class FlowCallback : public GRBCallback {
  public:
    double rootBound;

    FlowCallback(SVRPInstance &instance, Params &params, GRBModel &model,
                 ArcGRBVarMap &x, CVRPSEPSeparator *cvrpsepSeparator);
    void callback();

  private:
    double (GRBCallback::*solutionValue)(GRBVar);

    SVRPInstance &instance;
    Params &params;
    GRBModel &model;
    ArcGRBVarMap &x;
    CVRPSEPSeparator *cvrpsepSeparator;
    const int noImprovLimit = 5;
    const double noImprovEps = 1e-4;
    double prevDualBound;
    int prevNodeCount;
    int countNoImprov;

    void updateImprovData(int nodeCount, double dualBound);
};

#endif // FLOWCALLBACK_H
