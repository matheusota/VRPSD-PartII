#include "basiccallback.h"

BasicCallback::BasicCallback(SVRPInstance &instance, Params &params,
                             GRBModel &model, EdgeGRBVarMap &x,
                             std::vector<GRBVar> &scenariosVars,
                             CVRPSEPSeparator *cvrpsepSeparator,
                             GendreauSeparator *gendreauSeparator)
    : instance(instance), params(params), model(model), x(x),
      scenariosVars(scenariosVars) {
    this->cvrpsepSeparator = cvrpsepSeparator;
    this->gendreauSeparator = gendreauSeparator;
    prevNodeCount = -1;
    prevDualBound = 0.0;
    countNoImprov = 0;
    rootBound = 1e9;
}

void BasicCallback::updateImprovData(int nodeCount, double dualBound) {
    // New node, reset all data.
    if (nodeCount != prevNodeCount) {
        prevNodeCount = nodeCount;
        prevDualBound = dualBound;
        countNoImprov = 0;
        return;
    }

    // Update no improvements counter.
    if (dualBound - prevDualBound >= noImprovEps) {
        countNoImprov = 0;
    } else {
        countNoImprov++;
    }

    prevDualBound = dualBound;
}

void BasicCallback::callback() {
    try {
        int currNodeCount;
        double currSolCost = 0.0;
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
            solutionValue = &BasicCallback::getSolution;
            currSolCost = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
            currDualBound = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
        }
        // Fractional solution.
        else {
            assert(isMIPNode);
            currNodeCount =
                static_cast<int>(getDoubleInfo(GRB_CB_MIPNODE_NODCNT));
            solutionValue = &BasicCallback::getNodeRel;
            currDualBound = getDoubleInfo(GRB_CB_MIPNODE_OBJBND);
        }

        if (currNodeCount == 0) {
            rootBound = currDualBound;
        }
        updateImprovData(currNodeCount, currDualBound);

        // Construct xValue.
        EdgeValueMap xValue(instance.g, 0.0);
        for (EdgeIt e(instance.g); e != INVALID; ++e) {
            xValue[e] = (this->*solutionValue)(x[e]);
        }

        // Construct recourseValue.
        double totalRecourse = 0.0;
        std::vector<double> scenariosRecourse(instance.nScenarios, 0.0);
        for (int i = 0; i < instance.nScenarios; i++) {
            scenariosRecourse[i] = (this->*solutionValue)(scenariosVars[i]);
            totalRecourse += scenariosRecourse[i];
        }

        // Call separation routines.
        std::vector<CutData> separatedCuts;
        cvrpsepSeparator->separateCVRPSEPCuts(xValue, separatedCuts);

        if (gendreauSeparator != nullptr && separatedCuts.empty()) {
            gendreauSeparator->separateAggregatedGendreauCut(
                xValue, totalRecourse, separatedCuts);
        }

        // Add cuts.
        for (CutData &cutData : separatedCuts) {
            double sign = (cutData.sense == GRB_LESS_EQUAL) ? 1 : -1;
            if (sign * (cutData.LHS - cutData.RHS) <= 1e-6) {
                continue;
            }

            GRBLinExpr expr = 0.0;
            for (auto &[e, coef] : cutData.edgePairs) {
                expr += coef * x[e];
                assert(coef > -INFINITY && coef < INFINITY);
            }

            for (auto &[scenarioId, coef] : cutData.scenarioPairs) {
                expr += coef * scenariosVars[scenarioId];
            }

            if (cutData.name.find("GendreauCut") != std::string::npos) {
                for (int s = 0; s < instance.nScenarios; s++) {
                    expr += scenariosVars[s];
                }
            }

            assert(cutData.RHS > -INFINITY && cutData.RHS < INFINITY);

            // std::cout << cutData.name << std::endl;
            // std::cout << expr << " " << cutData.sense << " " << cutData.RHS
            //           << std::endl;
            addLazy(expr, cutData.sense, cutData.RHS);
            // lazyVector.push_back(std::make_tuple(expr, cutData.sense,
            //                                      cutData.RHS,
            //                                      cutData.name));
        }
    } catch (GRBException e) {
        std::cout << "BasicCallback: Error when calling callback. Error code: "
                  << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        exit(1);
    } catch (...) {
        std::cout << "Unknown error!" << std::endl;
        exit(1);
    }
}