#ifndef PARTIALROUTE_H
#define PARTIALROUTE_H

#include "mygraphlib.h"

class SVRPInstance;
struct PartialRouteEntry {
    std::vector<int> vertices;
    double recourseCost;
    double xInside = 0.0;
    double xLeft = 0.0;
    double xRight = 0.0;
    std::vector<double> scenarioDemands;

    PartialRouteEntry(std::vector<int> vertices, double recourseCost,
                      double xInside, double xLeft, double xRight,
                      int nScenarios)
        : vertices(vertices), recourseCost(recourseCost), xInside(xInside),
          xLeft(xLeft), xRight(xRight) {
        scenarioDemands.reserve(nScenarios);
    };
};

class PartialRoute {
  public:
    std::vector<PartialRouteEntry> entries;
    bool isRoute = false;
    double excess = 0.0; // 1 + (x(H) - |H| + 1)
    std::vector<double> totalScenarioDemands;

    PartialRoute(const SVRPInstance &instance);
    void buildFromRoute(const std::vector<int> &route);

  private:
    const SVRPInstance &instance;
};

#endif // PARTIALROUTE_H