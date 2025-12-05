#ifndef NODESCENARIOCALLBACK_H
#define NODESCENARIOCALLBACK_H

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
#include "sriseparator.h"
#include "cutdata.h"

class NodeScenariosCallback : public GRBCallback {
  public:
    double rootBound;

    NodeScenariosCallback(const SVRPInstance &instance, const Params &params,
                          const EdgeGRBVarMap &x, const NodeVectorGRBVarMap &y,
                          CVRPSEPSeparator &cvrpsepSeparator,
                          SRISeparator &sriSeparator);
    void callback();

  private:
    double (GRBCallback::*solutionValue)(GRBVar);

    const SVRPInstance &instance;
    const Params &params;
    const EdgeGRBVarMap &x;
    const NodeVectorGRBVarMap &y;
    CVRPSEPSeparator &cvrpsepSeparator;
    SRISeparator &sriSeparator;
};

#endif // NODESCENARIOCALLBACK_H
