#ifndef PartialRoutesBuilder_H
#define PartialRoutesBuilder_H

#include "mygraphlib.h"
#include "svrpinstance.h"
#include "graphviewer.h"
#include <cassert>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/dfs.h>
#include <list>
#include <utility>
#include <set>
#include <algorithm>
#include <random>
#include <chrono>
#include "params.h"
#include "partialroute.h"

typedef Dfs<Graph> DFS;

struct hash_vector {
    std::size_t operator()(const std::vector<int> &customers) const {
        std::vector<int> tmp = customers;
        std::sort(tmp.begin(), tmp.end());
        return boost::hash_range(tmp.begin(), tmp.end());
    }
};

class PartialRoutesBuilder {
  public:
    PartialRoutesBuilder(const SVRPInstance &instance, const Params &params,
                         const EdgeValueMap &xValue);

    const std::vector<PartialRoute> &getPartialRoutes() {
        return partialRoutes;
    }
    const std::unordered_set<std::vector<int>, hash_vector> &getCustomerSets() {
        return customerSets;
    }

  private:
    std::vector<PartialRoute> partialRoutes;
    std::unordered_set<std::vector<int>, hash_vector> customerSets;

    const SVRPInstance &instance;
    const Params &params;
    const EdgeValueMap &xValue;
    double EpsForIntegrality = 1e-4;

    Graph supportGraph;
    EdgeValueMap supportXValue;
    NodeValueMap edgeDepotValue;

    Graph tree; // this is actually a forest!
    NodeVectorIntMap treeVerticesMap;
    NodeBoolMap treeCutVertexMap;
    NodeValueMap treeRecourseCostMap;
    NodeValueMap treeNodeXValue; // x(S_i)
    NodeValueMap treeDemand;

    void buildPartialRoutes();
    void buildSupportGraph();
    void buildBlockCutTree();
    void pushNodeIntoPartialRoute(PartialRoute &partialRoute, Node treeNode);
    bool isRoute(const PartialRoute &partialRoute);
    void checkSortFillPartialRoute(PartialRoute &partialRoute);
};

#endif // PARTIALROUTESBUILDER_H
