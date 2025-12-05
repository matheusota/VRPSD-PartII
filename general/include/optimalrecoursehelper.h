#ifndef OPTIMALRECOURSEHELPER_H
#define OPTIMALRECOURSEHELPER_H

#include "mygraphlib.h"
#include <memory>
#include <vector>
#include <cassert>
#include <utility>
#include <queue>
#include <boost/functional/hash.hpp>
#include "unordered_set"
#include "partialroute.h"

class SVRPInstance;

class OptimalRecourseHelper {
  public:
    OptimalRecourseHelper(const SVRPInstance &instance_);
    ~OptimalRecourseHelper();

    double getRecourseLowerBound(const std::vector<int> &customers,
                                 int nVehicles) const;
    double getRecourseCost(const std::vector<std::vector<int>> &routes) const;
    double getRouteRecourseCost(const std::vector<int> &route) const;
    double getRouteRecourseCostInScenarioWithLP(
        const std::vector<int> &route, int scenarioId,
        std::vector<std::unordered_set<int>> &tightSubroutes) const;
    double getRouteRecourseCostInScenarioWithDP(const std::vector<int> &route,
                                                int scenarioId) const;
    double getPartialRouteRecourseCost(const PartialRoute &partialRoute) const;
    double getPartialRouteRecourseCostInScenarioWithDP(
        const PartialRoute &partialRoute, int scenarioId) const;

  private:
    const SVRPInstance &instance;
};
#endif // OPTIMALRECOURSEHELPER_H
