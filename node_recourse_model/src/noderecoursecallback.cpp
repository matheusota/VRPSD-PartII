#include "noderecoursecallback.h"

NodeRecourseCallback::NodeRecourseCallback(
    const SVRPInstance &instance, const Params &params, const EdgeGRBVarMap &x,
    const NodeGRBVarMap &y, const NodeRecourseCutBuilder &cutBuilder,
    CVRPSEPSeparator &cvrpsepSeparator, GendreauSeparator &gendreauSeparator,
    ParadaSeparator &paradaSeparator,
    PartialRouteSeparator &partialRouteSeparator,
    AggregatedSRISeparator &aggregatedSRISeparator)
    : instance(instance), params(params), x(x), y(y), cutBuilder(cutBuilder),
      cvrpsepSeparator(cvrpsepSeparator), gendreauSeparator(gendreauSeparator),
      paradaSeparator(paradaSeparator),
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
    // std::cout << "before checking" << std::endl;
    // cvrpsepSeparator.checkMyOldCuts();
    // std::cout << "after checking" << std::endl;
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
            assert(x[e].get(GRB_StringAttr_VarName).at(0) == 'x');
            xValue[e] = (this->*solutionValue)(x[e]);
        }

        // Construct yValue.
        NodeValueMap yValue(instance.g, 0.0);
        for (NodeIt v(instance.g); v != INVALID; ++v) {
            if (instance.g.id(v) != instance.depot) {
                assert(y[v].get(GRB_StringAttr_VarName).at(0) == 'y');
                yValue[v] = (this->*solutionValue)(y[v]);
            }
        }

        // Call separation routines.
        std::vector<CutData> separatedCuts;
        separateCuts(xValue, yValue, separatedCuts, isMIPSol, currNodeCount);

        // std::vector<size_t> cutOrder(separatedCuts.size());
        // std::iota(cutOrder.begin(), cutOrder.end(), 0);
        // std::sort(cutOrder.begin(), cutOrder.end(), [&](size_t i, size_t j) {
        //     return separatedCuts[i].violation > separatedCuts[j].violation;
        // });
        // int limit = std::min(static_cast<int>(separatedCuts.size()), 5);

        // Add cuts.
        for (const CutData &cutData : separatedCuts) {
            // CutData &cutData = separatedCuts[cutOrder[i]];
            assert(std::isfinite(cutData.LHS));
            assert(std::isfinite(cutData.RHS));
            double sign = (cutData.sense == GRB_LESS_EQUAL) ? 1 : -1;
            if (sign * (cutData.LHS - cutData.RHS) <= 1e-6) {
                continue;
            }

            addLazy(cutBuilder.buildCutExpr(cutData), cutData.sense,
                    cutData.RHS);
            // if (cutData.name == "CVRPSEP" ||
            //     cutData.name == "PartialRouteCut") {
            //     addLazy(cutBuilder.buildCutExpr(cutData), cutData.sense,
            //             cutData.RHS);
            // } else if (isMIPNode) {
            //     addCut(cutBuilder.buildCutExpr(cutData), cutData.sense,
            //            cutData.RHS);
            // }
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
                paradaSeparator.separateSetCut(xValue, yValue, customers,
                                               separatedCuts);
            }
        }
    }

    if (!separatedCuts.empty()) {
        return;
    }

    PartialRoutesBuilder partialRoutesBuilder(instance, params, xValue, yValue);

    if (params.paradaSet || params.sriCuts) {
        // Note these customer sets already contain one set of each connected
        // component.
        const auto &prCustomerSets = partialRoutesBuilder.getCustomerSets();
        for (const auto &customers : prCustomerSets) {
            if (params.sriCuts) {
                aggregatedSRISeparator.addCutFromSet(xValue, yValue, customers,
                                                     separatedCuts);
            } else if (params.paradaSet) {
                paradaSeparator.separateSetCut(xValue, yValue, customers,
                                               separatedCuts);
            }
        }
    }

    // Construct partial routes and try to separate more cuts.
    if (params.partialRouteCuts) {
        partialRouteSeparator.separatePartialRouteCuts(
            xValue, yValue, partialRoutesBuilder.getPartialRoutes(),
            separatedCuts);
    }

    if (!separatedCuts.empty()) {
        return;
    }

    if (!isInteger && params.sriFlowSeparation && params.sriCuts &&
        nodeCount <= 1e6) {
        aggregatedSRISeparator.flowHeuristicSeparation(xValue, yValue,
                                                       separatedCuts);
    }
    if (!separatedCuts.empty()) {
        return;
    }

    // Separate Gendreau cuts.
    if (params.gendreauCuts && separatedCuts.empty()) {
        double totalRecourse = 0.0;
        for (NodeIt v(instance.g); v != INVALID; ++v) {
            if (instance.g.id(v) != instance.depot) {
                totalRecourse += yValue[v];
            }
        }

        gendreauSeparator.separateAggregatedGendreauCut(xValue, totalRecourse,
                                                        separatedCuts);
    }
}