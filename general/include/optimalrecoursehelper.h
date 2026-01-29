#ifndef OPTIMALRECOURSEHELPER_H
#define OPTIMALRECOURSEHELPER_H

#include "mygraphlib.h"
#include <memory>
#include <vector>
#include <cassert>
#include <utility>
#include <queue>
#include <boost/functional/hash.hpp>
#include <chrono>
#include "unordered_set"
#include "partialroute.h"

class SVRPInstance;

class OptimalRecourseHelper {
  public:
    OptimalRecourseHelper(const SVRPInstance &instance_);
    ~OptimalRecourseHelper();
    mutable double totalDPTime = 0.0;
    mutable double totalLPTime = 0.0;

    double getRecourseLowerBound(const std::vector<int> &customers,
                                 int nVehicles) const;
    double getRecourseCost(const std::vector<std::vector<int>> &routes) const;
    double getPartialRouteRecourseCostWithLP(
        const PartialRoute &partialRoute,
        std::unordered_map<std::pair<int, int>, double,
                           boost::hash<std::pair<int, int>>> &betaDuals,
        std::unordered_map<std::tuple<int, int, int>, double,
                           boost::hash<std::tuple<int, int, int>>> &alphaDuals)
        const;
    double getPartialRouteRecourseCost(const PartialRoute &partialRoute) const;
    double getPartialRouteRecourseCost(
        const PartialRoute &partialRoute,
        std::unordered_map<std::pair<int, int>, double,
                           boost::hash<std::pair<int, int>>> &betaDuals,
        std::unordered_map<std::tuple<int, int, int>, double,
                           boost::hash<std::tuple<int, int, int>>> &alphaDuals)
        const;
    double getPartialRouteRecourseCostInScenarioWithDP(
        const PartialRoute &partialRoute, int scenarioId) const;

  private:
    const SVRPInstance &instance;
    GRBEnv env;
};
#endif // OPTIMALRECOURSEHELPER_H
