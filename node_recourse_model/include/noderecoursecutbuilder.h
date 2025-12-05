#ifndef NODERECOURSECUTBUILDER_H
#define NODERECOURSECUTBUILDER_H

#include <gurobi_c++.h>
#include "mygraphlib.h"
#include "svrpinstance.h"
#include "params.h"
#include <cassert>
#include <lemon/list_graph.h>
#include <list>
#include <utility>
#include <ctime> // For CLOCKS_PER_SEC
#include <chrono>
#include <cmath>
#include "cutdata.h"

class NodeRecourseCutBuilder {
  public:
    NodeRecourseCutBuilder(const SVRPInstance &instance, const Params &params,
                           const EdgeGRBVarMap &x, const NodeGRBVarMap &y);
    const GRBLinExpr buildCutExpr(const CutData &cutData) const;

  private:
    const SVRPInstance &instance;
    const Params &params;
    const EdgeGRBVarMap &x;
    const NodeGRBVarMap &y;
};

#endif // NODERECOURSECUTBUILDER_H
