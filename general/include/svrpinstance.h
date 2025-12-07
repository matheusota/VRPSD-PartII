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
    int nScenarios;
    std::vector<std::vector<int>> scenariosMatrix;
    int k, n, m, depot;
    std::vector<std::pair<double, int>> sortedScenarios;
    double capacity;
    ClassicalRecourseHelper classicalRecourseHelper;
    OptimalRecourseHelper optimalRecourseHelper;

    SVRPInstance();
    void print() const;
    void sortScenarios();
    bool checkDeterministicRoute(const std::vector<int> &route) const;
    double getEdgeRecourseCost(Node v) const;
};
#endif // SVRPINSTANCE_H
