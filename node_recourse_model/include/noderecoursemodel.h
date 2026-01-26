#ifndef NODERECOURSEMODEL_H
#define NODERECOURSEMODEL_H

#include <gurobi_c++.h>
#include <utility>
#include <map>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <stack>
#include <chrono>
#include <cmath>
#include "mygraphlib.h"
#include "svrpinstance.h"
#include "graphviewer.h"
#include "params.h"
#include "cvrpsepseparator.h"
#include "paradaseparator.h"
#include "noderecoursecallback.h"
#include "partialroute.h"
#include "partialroutesbuilder.h"
#include "partialrouteseparator.h"
#include "aggregatedsriseparator.h"
#include "cutpool.h"
#include "noderecoursecutbuilder.h"
#include <boost/functional/hash.hpp>

class NodeRecourseModel {
  public:
    NodeRecourseModel(SVRPInstance &instance, const Params &params,
                      CutPool *cutPool_ = nullptr);

    bool solve(SVRPSolution &solution);
    void setWarmStart(const SVRPSolution &primalSolution);

  private:
    const SVRPInstance &instance;
    const Params &params;
    CutPool *cutPool;
    double initialRootBound;

    // Gurobi stuff.
    GRBEnv env;
    GRBModel model;
    EdgeGRBVarMap x;
    NodeGRBVarMap y;
    NodeRecourseCutBuilder cutBuilder;

    CVRPSEPSeparator cvrpsepSeparator;
    ParadaSeparator paradaSeparator;
    PartialRouteSeparator partialRouteSeparator;
    AggregatedSRISeparator aggregatedSRISeparator;
    NodeRecourseCallback callback;

    void setBasicModel();
    void solveRootLP(
        std::chrono::time_point<std::chrono::high_resolution_clock> &started);
    void getCurrentSolution(EdgeValueMap &xValue, NodeValueMap &yValue);
    bool addSeparatedCuts(const std::vector<CutData> &separatedCuts);
    void setSolution(SVRPSolution &solution);
    void addInitialCuts(const std::vector<CutData> &initialCuts);
    void cleanupCuts();
};

#endif // NODERECOURSEMODEL_H