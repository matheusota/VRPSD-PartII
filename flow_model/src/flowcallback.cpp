#include "flowcallback.h"

FlowCallback::FlowCallback(SVRPInstance &instance, Params &params,
                           GRBModel &model, ArcGRBVarMap &x,
                           CVRPSEPSeparator *cvrpsepSeparator)
    : instance(instance), params(params), model(model), x(x) {
    this->cvrpsepSeparator = cvrpsepSeparator;
    prevNodeCount = -1;
    prevDualBound = 0.0;
    countNoImprov = 0;
    rootBound = 1e9;
}

void FlowCallback::updateImprovData(int nodeCount, double dualBound) {
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

void FlowCallback::callback() {
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
            currNodeCount = (int)getDoubleInfo(GRB_CB_MIPSOL_NODCNT);
            solutionValue = &FlowCallback::getSolution;
            currSolCost = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
            currDualBound = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
        }
        // Fractional solution.
        else {
            assert(isMIPNode);
            currNodeCount = (int)getDoubleInfo(GRB_CB_MIPNODE_NODCNT);
            solutionValue = &FlowCallback::getNodeRel;
            currDualBound = getDoubleInfo(GRB_CB_MIPNODE_OBJBND);
        }

        if (currNodeCount == 0) {
            rootBound = currDualBound;
        }
        updateImprovData(currNodeCount, currDualBound);

        // Construct xValue.
        ArcValueMap xValue(instance.d_g, 0.0);
        for (ArcIt a(instance.d_g); a != INVALID; ++a) {
            xValue[a] = (this->*solutionValue)(x[a]);
        }

        // Call separation routines.
        std::vector<CutData> separatedCuts;
        cvrpsepSeparator->separateCVRPSEPCuts(xValue, separatedCuts);

        // Add cuts.
        for (CutData &cutData : separatedCuts) {
            double sign = (cutData.sense == GRB_LESS_EQUAL) ? 1 : -1;
            if (sign * (cutData.LHS - cutData.RHS) <= 1e-6) {
                continue;
            }

            GRBLinExpr expr = 0.0;
            for (std::pair<Edge, double> &p : cutData.edgePairs) {
                std::pair<Arc, Arc> arcPair = instance.getArcsFromEdge(p.first);

                expr +=
                    p.second * x[arcPair.first] + p.second * x[arcPair.second];
                assert(p.second > -INFINITY && p.second < INFINITY);
            }

            assert(cutData.RHS > -INFINITY && cutData.RHS < INFINITY);

            addLazy(expr, cutData.sense, cutData.RHS);
            // lazyVector.push_back(std::make_tuple(expr, cutData.sense,
            //                                      cutData.RHS, cutData.name));
        }
    } catch (GRBException e) {
        std::cout << "FlowCallback: Error when calling callback. Error code: "
                  << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        exit(1);
    } catch (...) {
        std::cout << "Unknown error!" << std::endl;
        exit(1);
    }
}