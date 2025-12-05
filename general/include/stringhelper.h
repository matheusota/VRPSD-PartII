#ifndef STRINGHELPER_H
#define STRINGHELPER_H

#include "svrpinstance.h"
#include "svrpsolution.h"
#include <string>

namespace stringhelper {
std::string instanceAsString(const SVRPInstance &instance);
std::string instanceDescriptionAsString(const SVRPInstance &instance);
std::string nodesAndDemandAsString(const SVRPInstance &instance);
std::string edgesAndDemandAsString(const SVRPInstance &instance);
std::string vti(DNode v, const SVRPInstance &instance);
std::string demandAsString(const SVRPInstance &instance);
std::string solutionAsString(const SVRPInstance &instance,
                             SVRPSolution &solution);
std::string routesAsString(const SVRPSolution &solution);
std::string valuesAsString(const SVRPSolution &solution);
} // namespace stringhelper
#endif // STRINGHELPER_H
