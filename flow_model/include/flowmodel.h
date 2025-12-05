#ifndef FLOW_MODEL_H
#define FLOW_MODEL_H

#include <gurobi_c++.h>
#include <utility>
#include <map>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <stack>
#include "mygraphlib.h"
#include "svrpinstance.h"
#include "graphviewer.h"
#include "params.h"
#include "cvrpsepseparator.h"
#include "flowcallback.h"

class FlowModel {
  public:
    FlowModel(SVRPInstance &instance, Params &params);
    ~FlowModel();
    bool solve(SVRPSolution &solution);

  private:
    SVRPInstance &instance;
    Params &params;

    // Gurobi stuff.
    GRBEnv env;
    GRBModel model;
    ArcGRBVarMap x;
    std::vector<DNodeGRBVarMap *> y;
    std::vector<DNodeGRBVarMap *> z;
    std::vector<ArcGRBVarMap *> f;
    CVRPSEPSeparator *cvrpsepSeparator;
    FlowCallback *callback;

    void setBasicModel();
    void solveRootLP();
    bool addSeparatedCuts(std::vector<CutData> &separatedCuts);
    void setSolution(SVRPSolution &solution);
};

#endif // FLOW_MODEL_H