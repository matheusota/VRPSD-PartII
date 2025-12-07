#ifndef NODERECOURSECALLBACK_H
#define NODERECOURSECALLBACK_H

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
#include <cmath>
#include "cvrpsepseparator.h"
#include "paradaseparator.h"
#include "partialroute.h"
#include "partialroutesbuilder.h"
#include "partialrouteseparator.h"
#include "aggregatedsriseparator.h"
#include "cutdata.h"
#include "noderecoursecutbuilder.h"

class NodeRecourseCallback : public GRBCallback {
  public:
    double rootBound;

    NodeRecourseCallback(const SVRPInstance &instance, const Params &params,
                         const EdgeGRBVarMap &x, const NodeGRBVarMap &y,
                         const NodeRecourseCutBuilder &cutBuilder,
                         CVRPSEPSeparator &cvrpsepSeparator,
                         ParadaSeparator &paradaSeparator,
                         PartialRouteSeparator &partialRouteSeparator,
                         AggregatedSRISeparator &aggregatedSRISeparator);
    void callback();

    void separateCuts(const EdgeValueMap &xValue, const NodeValueMap &yValue,
                      std::vector<CutData> &separatedCuts, bool isInteger,
                      int nodeCount);

  private:
    double (GRBCallback::*solutionValue)(GRBVar);

    const SVRPInstance &instance;
    const Params &params;
    const EdgeGRBVarMap &x;
    const NodeGRBVarMap &y;
    const NodeRecourseCutBuilder &cutBuilder;
    CVRPSEPSeparator &cvrpsepSeparator;
    ParadaSeparator &paradaSeparator;
    PartialRouteSeparator &partialRouteSeparator;
    AggregatedSRISeparator &aggregatedSRISeparator;
};

#endif // NODERECOURSECALLBACK_H
