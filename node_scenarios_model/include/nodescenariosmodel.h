#ifndef NODESCENARIOMODEL_H
#define NODESCENARIOMODEL_H

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
#include <algorithm>
#include <numeric>
#include "mygraphlib.h"
#include "svrpinstance.h"
#include "graphviewer.h"
#include "params.h"
#include "cvrpsepseparator.h"
#include "sriseparator.h"
#include "cutpool.h"

class NodeScenariosModel {
  public:
    NodeScenariosModel(const SVRPInstance &instance, const Params &params,
                       CutPool &cutPool);

    bool solve(SVRPSolution &solution);
    void setPrimalSolution(const SVRPSolution &primalSolution);

  private:
    const SVRPInstance &instance;
    const Params &params;
    CutPool &cutPool;

    // Gurobi stuff.
    GRBEnv env;
    GRBModel model;
    EdgeGRBVarMap x;
    NodeVectorGRBVarMap y;

    // Parameters for In/Out separation.
    EdgeValueMap interiorX;
    NodeVectorValueMap interiorY;
    double LB = 0.0;
    double UB = 1e9;
    double lambda = 1.0;

    CVRPSEPSeparator cvrpsepSeparator;
    SRISeparator sriSeparator;

    void setBasicModel();
    void solveRootLP(
        std::chrono::time_point<std::chrono::high_resolution_clock> &started);
    int addSeparatedCuts(EdgeValueMap &xValue, NodeVectorValueMap &yValue,
                         EdgeValueMap &candX, NodeVectorValueMap &candY,
                         std::vector<CutData> &separatedCuts);
    double buildLagrangianCutAndGetDualObj();
    void solveLPAndUpdateSolutions(EdgeValueMap &xValue,
                                   NodeVectorValueMap &yValue,
                                   EdgeValueMap &candX,
                                   NodeVectorValueMap &candY);
    void cleanupCuts();
};

#endif // NODESCENARIOMODEL_H