#ifndef GRAPHVIEWER_H
#define GRAPHVIEWER_H
#include <iomanip> // setprecision
#include <sstream> // stringstream
#include "mygraphlib.h"
#include <string>
#include <vector>
#include "svrpinstance.h"
#include "svrpsolution.h"

namespace graphviewer {
void viewVRPSolution(const SVRPInstance &instance, const EdgeValueMap &x,
                     const std::string title);
void viewVRPSolution(const SVRPInstance &instance, const ArcValueMap &x,
                     const std::string title);
void viewVRPSolution(const SVRPInstance &instance, const SVRPSolution &solution,
                     const std::string title);
void viewSupportWithNodeAndEdgeMap(Graph &g, EdgeValueMap &x, NodePosMap &posx,
                                   NodePosMap &posy, NodeBoolMap &nodeMap,
                                   EdgeIntMap &edgeMap);
} // namespace graphviewer
#endif // GRAPHVIEWER_H
