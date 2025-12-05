#ifndef CLASSICALRECOURSEHELPER_H
#define CLASSICALRECOURSEHELPER_H

#include "mygraphlib.h"
#include <memory>
#include <vector>
#include <cassert>
#include <utility>
#include <queue>
#include "partialroute.h"

class SVRPInstance;

class ClassicalRecourseHelper {
  public:
    ClassicalRecourseHelper(const SVRPInstance &instance_);
    ~ClassicalRecourseHelper();

    double getRecourseCost(const std::vector<std::vector<int>> &routes) const;
    double getRouteRecourseCost(const std::vector<int> &route) const;
    double getRouteRecourseCostInScenario(const std::vector<int> &route,
                                          int scenarioId, bool reverse) const;
    double
    getRecourseCostInScenario(const std::vector<std::vector<int>> &routes,
                              int scenarioId) const;
    bool isRouteReversed(const std::vector<int> &route) const;

    double
    getPartialRouteRecourseAndFailureNodes(const PartialRoute &partialRoute,
                                           NodeValueMap &nodeRecourseMap) const;

    double getRouteRecourseCostAndFailures(const std::vector<int> &route,
                                           NodeVectorValueMap &failures) const;

  private:
    const SVRPInstance &instance;
};
#endif // OPTIMALRECOURSEHELPER_H
