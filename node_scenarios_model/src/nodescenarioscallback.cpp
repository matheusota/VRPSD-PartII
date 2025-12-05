#include "nodescenarioscallback.h"

NodeScenariosCallback::NodeScenariosCallback(const SVRPInstance &instance,
                                             const Params &params,
                                             const EdgeGRBVarMap &x,
                                             const NodeVectorGRBVarMap &y,
                                             CVRPSEPSeparator &cvrpsepSeparator,
                                             SRISeparator &sriSeparator)
    : instance(instance), params(params), x(x), y(y),
      cvrpsepSeparator(cvrpsepSeparator), sriSeparator(sriSeparator) {
    rootBound = 1e9;
}

void NodeScenariosCallback::callback() {
    try {
        int currNodeCount;
        double currDualBound = 0.0;
        bool isMIPSol = (where == GRB_CB_MIPSOL);
        bool isMIPNode = (where == GRB_CB_MIPNODE) &&
                         (getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL);

        if (!isMIPSol && !isMIPNode) {
            return;
        }

        // Get solution values according to solution type.
        // Integer solution.
        if (isMIPSol) {
            currNodeCount =
                static_cast<int>(getDoubleInfo(GRB_CB_MIPSOL_NODCNT));
            solutionValue = &NodeScenariosCallback::getSolution;
            currDualBound = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
        }
        // Fractional solution.
        else {
            assert(isMIPNode);
            currNodeCount =
                static_cast<int>(getDoubleInfo(GRB_CB_MIPNODE_NODCNT));
            solutionValue = &NodeScenariosCallback::getNodeRel;
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

        // Construct yValue and totalRecourse.
        double totalRecourse = 0.0;
        NodeVectorValueMap yValue(
            instance.g, std::vector<double>(instance.nScenarios, 0.0));
        NodeValueMap recourseValue(instance.g, 0.0);
        for (int scenarioId = 0; scenarioId < instance.nScenarios;
             scenarioId++) {
            for (NodeIt v(instance.g); v != INVALID; ++v) {
                if (instance.g.id(v) != instance.depot) {
                    yValue[v][scenarioId] =
                        (this->*solutionValue)(y[v][scenarioId]);
                    recourseValue[v] =
                        instance.getEdgeRecourseCost(v) * yValue[v][scenarioId];
                    totalRecourse += recourseValue[v];
                }
            }
        }

        // Call separation routines.
        std::vector<CutData> separatedCuts;
        cvrpsepSeparator.separateCVRPSEPCuts(xValue, separatedCuts);

        if (params.sriCuts && separatedCuts.empty() && isMIPSol) {
            sriSeparator.integerSeparation(xValue, yValue, separatedCuts);
            // int addedCuts = nodeRecourseSeparator->heuristicSeparation(
            //     xValue, yValue, separatedCuts);

            // if (addedCuts == 0 && isMIPSol) {
            //     nodeRecourseSeparator->integerSeparation(xValue, yValue,
            //                                              separatedCuts);
            // }

            // if (addedCuts == 0 && params.mipSeparation && isMIPSol) {
            //     nodeRecourseSeparator->mipSeparation(xValue, yValue,
            //                                          separatedCuts);
            // }
        }

        // Add cuts.
        for (CutData &cutData : separatedCuts) {
            double sign = (cutData.sense == GRB_LESS_EQUAL) ? 1 : -1;
            if (sign * (cutData.LHS - cutData.RHS) <= 1e-6) {
                continue;
            }

            GRBLinExpr expr = 0.0;
            for (auto [e, coef] : cutData.edgePairs) {
                assert(coef > -INFINITY && coef < INFINITY);
                expr += coef * x[e];
            }

            for (auto [scenarioId, scenarioCoef] : cutData.scenarioPairs) {
                for (auto [nodeId, nodeCoef] : cutData.nodePairs) {
                    assert(nodeId != instance.depot);
                    Node v = instance.g.nodeFromId(nodeId);
                    expr += nodeCoef * y[v][scenarioId];
                }
            }

            assert(cutData.RHS > -INFINITY && cutData.RHS < INFINITY);

            addLazy(expr, cutData.sense, cutData.RHS);
        }
    } catch (GRBException e) {
        std::cout << "NodeScenariosCallback: Error when calling callback. "
                     "Error code: "
                  << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        exit(1);
    } catch (...) {
        std::cout << "Unknown error!" << std::endl;
        exit(1);
    }
}