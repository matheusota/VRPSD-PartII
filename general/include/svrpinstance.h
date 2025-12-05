#ifndef SVRPINSTANCE_H
#define SVRPINSTANCE_H

#include "classicalrecoursehelper.h"
#include "optimalrecoursehelper.h"
#include "mygraphlib.h"
#include <memory>
#include <vector>
#include <cassert>
#include <utility>
#include <queue>

class SVRPInstance {
  public:
    Graph g;
    NodePosMap posx;
    NodePosMap posy;
    NodeValueMap demand;
    EdgeValueMap weight;
    NodeIntMap upperBoundMap;
    Digraph d_g;
    DNodePosMap d_posx;
    DNodePosMap d_posy;
    DNodeValueMap d_demand;
    ArcValueMap d_weight;
    int nScenarios;
    std::vector<std::vector<int>> scenariosMatrix;
    int k, n, m, nrows, depot;
    std::vector<std::pair<double, int>> sortedScenarios;
    double capacity;
    ClassicalRecourseHelper classicalRecourseHelper;
    OptimalRecourseHelper optimalRecourseHelper;

    SVRPInstance();
    void sortScenarios();
    Node getNodeFromDNode(DNode v) const;
    DNode getDNodeFromNode(Node v) const;
    Edge getEdgeFromArc(Arc a) const;
    std::pair<Arc, Arc> getArcsFromEdge(Edge e) const;
    bool checkDeterministicRoute(const std::vector<int> &route) const;
    double getEdgeRecourseCost(Node v) const;
    double getArcRecourseCost(DNode v) const;
};
#endif // SVRPINSTANCE_H
