#include "noderecoursecallback.h"

NodeRecourseCallback::NodeRecourseCallback(
    const SVRPInstance &instance, const Params &params, const EdgeGRBVarMap &x,
    const NodeGRBVarMap &y, const NodeRecourseCutBuilder &cutBuilder,
    CVRPSEPSeparator &cvrpsepSeparator, ParadaSeparator &paradaSeparator,
    PartialRouteSeparator &partialRouteSeparator,
    AggregatedSRISeparator &aggregatedSRISeparator)
    : instance(instance), params(params), x(x), y(y), cutBuilder(cutBuilder),
      cvrpsepSeparator(cvrpsepSeparator), paradaSeparator(paradaSeparator),
      partialRouteSeparator(partialRouteSeparator),
      aggregatedSRISeparator(aggregatedSRISeparator) {
    rootBound = 1e9;
}

void NodeRecourseCallback::callback() {
    int currNodeCount;
    double currDualBound = 0.0;
    bool isMIPSol = (where == GRB_CB_MIPSOL);
    bool isMIPNode = (where == GRB_CB_MIPNODE) &&
                     (getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL);
    if (!isMIPSol && !isMIPNode) {
        return;
    }

    try {
        // Get solution values according to solution type.
        // Integer solution.
        if (isMIPSol) {
            currNodeCount =
                static_cast<int>(getDoubleInfo(GRB_CB_MIPSOL_NODCNT));
            solutionValue = &NodeRecourseCallback::getSolution;
            currDualBound = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
        }
        // Fractional solution.
        else {
            assert(isMIPNode);
            currNodeCount =
                static_cast<int>(getDoubleInfo(GRB_CB_MIPNODE_NODCNT));
            solutionValue = &NodeRecourseCallback::getNodeRel;
            currDualBound = getDoubleInfo(GRB_CB_MIPNODE_OBJBND);
        }

        if (currNodeCount == 0) {
            rootBound = currDualBound;
        }

        // Construct xValue.
        EdgeValueMap xValue(instance.g, 0.0);
        for (EdgeIt e(instance.g); e != INVALID; ++e) {
            xValue[e] = (this->*solutionValue)(x[e]);
        }

        // Construct yValue.
        NodeValueMap yValue(instance.g, 0.0);
        for (NodeIt v(instance.g); v != INVALID; ++v) {
            if (instance.g.id(v) != instance.depot) {
                yValue[v] = (this->*solutionValue)(y[v]);
            }
        }

        // Call separation routines.
        std::vector<CutData> separatedCuts;
        separateCuts(xValue, yValue, separatedCuts, isMIPSol, currNodeCount);

        // Add cuts.
        for (const CutData &cutData : separatedCuts) {
            assert(std::isfinite(cutData.LHS));
            assert(std::isfinite(cutData.RHS));
            double sign = (cutData.sense == GRB_LESS_EQUAL) ? 1 : -1;
            if (sign * (cutData.LHS - cutData.RHS) <= 1e-6) {
                continue;
            }

            addLazy(cutBuilder.buildCutExpr(cutData), cutData.sense,
                    cutData.RHS);
        }
    } catch (GRBException e) {
        std::cout
            << "NodeRecourseCallback: Error when calling callback. Error code: "
            << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        exit(1);
    } catch (...) {
        std::cout << "Unknown error!" << std::endl;
        exit(1);
    }
}

void NodeRecourseCallback::separateCuts(const EdgeValueMap &xValue,
                                        const NodeValueMap &yValue,
                                        std::vector<CutData> &separatedCuts,
                                        bool isInteger, int nodeCount) {
    NodeValueMap recourseValue(instance.g, 0.0);
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) != instance.depot) {
            recourseValue[v] = (instance.getEdgeRecourseCost(v) /
                                static_cast<double>(instance.nScenarios)) *
                               yValue[v];
        }
    }

    // First separate CVRPSEP cuts.
    separatedCuts.clear();
    cvrpsepSeparator.separateCVRPSEPCuts(xValue, separatedCuts);
    const std::vector<std::vector<int>> &cvrpsepCustomerSets =
        cvrpsepSeparator.getCustomerSets();
    const std::vector<int> &cvrpsepNVehicles = cvrpsepSeparator.getNVehicles();

    // Try to separate cuts from CVRPSEP.
    if (params.sriCuts || params.paradaSet) {
        for (const auto &customers : cvrpsepCustomerSets) {
            if (params.sriCuts) {
                aggregatedSRISeparator.addCutFromSet(xValue, yValue, customers,
                                                     separatedCuts);
            } else if (params.paradaSet) {
                paradaSeparator.separateSetCut(xValue, recourseValue, customers,
                                               separatedCuts);
            }
        }
    }

    if (!separatedCuts.empty()) {
        return;
    }

    // Try to separate inequalities from partial routes.
    auto started = chrono::high_resolution_clock::now();
    PartialRoutesBuilder partialRoutesBuilder(instance, params, xValue);
    auto done = chrono::high_resolution_clock::now();
    partialRouteSeparator.time +=
        static_cast<double>(
            chrono::duration_cast<chrono::microseconds>(done - started)
                .count()) /
        1e6;

    if (params.paradaSet || params.sriCuts) {
        // Note these customer sets already contain one set for each connected
        // component.
        const auto &prCustomerSets = partialRoutesBuilder.getCustomerSets();
        for (const auto &customers : prCustomerSets) {
            if (params.sriCuts) {
                aggregatedSRISeparator.addCutFromSet(xValue, yValue, customers,
                                                     separatedCuts);
            } else if (params.paradaSet) {
                paradaSeparator.separateSetCut(xValue, recourseValue, customers,
                                               separatedCuts);
            }
        }
    }

    // Construct partial routes and try to separate more cuts.
    if (params.partialRouteCuts) {
        partialRouteSeparator.separatePartialRouteCuts(
            xValue, recourseValue, partialRoutesBuilder.getPartialRoutes(),
            separatedCuts);
    }

    if (!separatedCuts.empty()) {
        return;
    }

    // Lastly, we separate aggregated SRIs using the flow heuristic.
    if (!isInteger && params.sriFlowSeparation && params.sriCuts &&
        nodeCount <= 5000) {
        aggregatedSRISeparator.flowHeuristicSeparation(xValue, yValue,
                                                       separatedCuts);
    }
}