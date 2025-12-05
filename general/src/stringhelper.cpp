#include "stringhelper.h"

namespace stringhelper {
std::string instanceDescriptionAsString(const SVRPInstance &instance) {
    std::stringstream ss;

    ss << "n          : " << instance.n << endl;
    ss << "m          : " << instance.m << endl;
    ss << "k          : " << instance.k << endl;
    ss << "capacity   : " << instance.capacity << endl;
    ss << "depot      : " << instance.depot << endl;

    ss << "nodes     :" << endl;
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        ss << " " << instance.g.id(v) << ": " << instance.demand[v] << endl;
    }
    ss << endl;

    ss << "-------------------------" << endl;
    for (int j = 0; j < instance.nScenarios; j++) {
        ss << "SCENARIO " << j << " :";
        for (int i = 0; i < instance.n; i++) {
            ss << instance.scenariosMatrix[i][j] << " ";
        }
        ss << endl;
    }
    return ss.str();
}

std::string instanceAsString(const SVRPInstance &instance) {
    std::stringstream ss;

    ss << instanceDescriptionAsString(instance);
    ss << nodesAndDemandAsString(instance);
    ss << edgesAndDemandAsString(instance);
    ss << demandAsString(instance);

    return ss.str();
}

std::string nodesAndDemandAsString(const SVRPInstance &instance) {
    std::stringstream ss;
    ss << "nodes     :";
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        ss << " " << instance.g.id(v) << "." << instance.demand[v];
    }
    ss << endl;
    return ss.str();
}

std::string edgesAndDemandAsString(const SVRPInstance &instance) {
    std::stringstream ss;
    ss << "edges       :";
    for (EdgeIt a(instance.g); a != INVALID; ++a) {
        ss << " (" << instance.g.id(instance.g.u(a)) << ", "
           << instance.g.id(instance.g.v(a)) << "; " << instance.weight[a]
           << ") ";
    }
    ss << endl;
    return ss.str();
}

std::string demandAsString(const SVRPInstance &instance) {
    std::stringstream ss;
    ss << "demand      :";
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        ss << "[" << instance.g.id(v) << ": " << instance.demand[v] << "] ";
    }
    ss << endl;
    return ss.str();
}

std::string solutionAsString(const SVRPInstance &instance,
                             const SVRPSolution &sol) {
    std::stringstream ss;
    ss << instanceDescriptionAsString(instance);
    ss << valuesAsString(sol);
    ss << routesAsString(sol);
    return ss.str();
}

std::string routesAsString(const SVRPSolution &sol) {
    stringstream ss;

    ss << "Solution:" << endl;
    for (auto route : sol.routes) {
        for (int i : route) {
            ss << to_string(i) << " ";
        }
        ss << endl;
    }

    return ss.str();
}

std::string valuesAsString(const SVRPSolution &sol) {
    std::stringstream ss;
    ss << "lowerBound : " << sol.lowerBound << endl;
    ss << "cost       : " << sol.cost << endl;
    ss << "upperBound : " << sol.upperBound << endl;
    return ss.str();
}
} // namespace stringhelper
