#ifndef BASICMODEL_H
#define BASICMODEL_H

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
#include "mygraphlib.h"
#include "svrpinstance.h"
#include "graphviewer.h"
#include "params.h"
#include "cvrpsepseparator.h"
#include "gendreauseparator.h"
#include "basiccallback.h"

class BasicModel {
  public:
    BasicModel(SVRPInstance &instance, Params &params);
    ~BasicModel();

    bool solve(SVRPSolution &solution);

  private:
    SVRPInstance &instance;
    Params &params;

    // Gurobi stuff.
    GRBEnv env;
    GRBModel model;
    EdgeGRBVarMap x;
    std::vector<GRBVar> scenariosVars;
    BasicCallback *callback;

    CVRPSEPSeparator *cvrpsepSeparator;
    GendreauSeparator *gendreauSeparator;

    void setBasicModel();
    void solveRootLP(
        std::chrono::time_point<std::chrono::high_resolution_clock> &started);
    void getCurrentSolution(EdgeValueMap &xValue,
                            std::vector<double> &scenarioRecourses,
                            double &totalRecourse);
    bool addSeparatedCuts(std::vector<CutData> &separatedCuts);
    void setSolution(SVRPSolution &solution);
};

#endif // BASICMODEL_H