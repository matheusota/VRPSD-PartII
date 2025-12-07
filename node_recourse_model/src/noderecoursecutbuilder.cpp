#include "noderecoursecutbuilder.h"

NodeRecourseCutBuilder::NodeRecourseCutBuilder(const SVRPInstance &instance,
                                               const Params &params,
                                               const EdgeGRBVarMap &x,
                                               const NodeGRBVarMap &y)
    : instance(instance), params(params), x(x), y(y) {}

const GRBLinExpr
NodeRecourseCutBuilder::buildCutExpr(const CutData &cutData) const {
    GRBLinExpr expr = 0.0;

    for (const auto &[e, coef] : cutData.edgePairs) {
        assert(e != INVALID);
        assert(std::isfinite(coef));
        assert(x[e].get(GRB_StringAttr_VarName).at(0) == 'x');
        expr += coef * x[e];
    }

    for (const auto [id, coef] : cutData.nodePairs) {
        assert(id >= 1 && id < instance.n);
        assert(std::isfinite(coef));
        Node v = instance.g.nodeFromId(id);
        assert(y[v].get(GRB_StringAttr_VarName).at(0) == 'y');

        double newCoef = coef;
        if (cutData.name == "PartialRouteCut" || cutData.name == "SetCut") {
            newCoef *= instance.getEdgeRecourseCost(v) /
                       static_cast<double>(instance.nScenarios);
        }
        expr += newCoef * y[v];
    }

    return expr;
}