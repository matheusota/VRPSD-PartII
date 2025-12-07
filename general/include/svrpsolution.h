#ifndef SVRPSOLUTION_H
#define SVRPSOLUTION_H

#include "mygraphlib.h"
#include "svrpinstance.h"
#include <vector>

class SVRPSolution {
  public:
    std::vector<std::vector<int>> routes;
    double lowerBound = 0.0;
    double deterministicCost = 0.0;
    double recourseCost = 0.0;
    double optimalRecourseCost = 0.0;
    double modelRecourseCost = 0.0;
    double cost = 0.0;
    double upperBound = 1e9;
    double rootBound = 0.0;
    double time = 0.0;
    int cvrpsepCuts = 0;
    double cvrpsepTime = 0.0;
    int paradaPathCuts = 0;
    int paradaSetCuts = 0;
    double paradaTime = 0.0;
    int partialRouteCuts = 0;
    double partialRouteTime = 0.0;
    int sriCuts = 0;
    double sriTime = 0.0;
    int aggregatedSriCuts = 0;
    double aggregatedSriTime = 0.0;
    int nodesExplored = 0;
    double gap = 0.0;
    double rootTime = 0.0;
    bool solved = false;

    void print() const;
    void save(std::string outputName) const;
    void setRoutesFromSolution(const SVRPInstance &instance,
                               EdgeIntMap &edgeCount);
    void setRoutesFromSolution(const SVRPInstance &instance,
                               ArcIntMap &edgeCount);
    std::vector<Edge> getEdgesInSolution(const SVRPInstance &instance) const;
};
#endif // SVRPSOLUTION_H
