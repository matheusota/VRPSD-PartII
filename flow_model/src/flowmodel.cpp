#include "flowmodel.h"

FlowModel::FlowModel(SVRPInstance &instance, Params &params)
    : instance(instance), params(params), model(env), x(instance.d_g) {
    cvrpsepSeparator = new CVRPSEPSeparator(instance, params);

    // Set Gurobi callback.
    callback = new FlowCallback(instance, params, model, x, cvrpsepSeparator);
    model.setCallback(callback);

    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        y.push_back(new DNodeGRBVarMap(instance.d_g));
        z.push_back(new DNodeGRBVarMap(instance.d_g));
        f.push_back(new ArcGRBVarMap(instance.d_g));
    }

    setBasicModel();
}

FlowModel::~FlowModel() {
    delete cvrpsepSeparator;
    delete callback;

    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        delete z[scenarioId];
        delete f[scenarioId];
    }
}

void FlowModel::setBasicModel() {
    model.set(GRB_StringAttr_ModelName, "Flow2Model");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Must set to use lazy contraints.
    model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
    model.getEnv().set(GRB_IntParam_PreCrush, 1);

    // Some other parameters.
    model.getEnv().set(GRB_DoubleParam_TimeLimit, params.timeLimit);
    model.getEnv().set(GRB_IntParam_Threads, 4);
    model.getEnv().set(GRB_DoubleParam_MIPGap, 1e-6);
    model.getEnv().set(GRB_IntParam_Cuts, 3);
    model.getEnv().set(GRB_IntParam_Presolve, 2);

    // Declare x variables.
    for (ArcIt a(instance.d_g); a != INVALID; ++a) {
        x[a] = model.addVar(
            0.0, 1.0, instance.d_weight[a], GRB_CONTINUOUS,
            "x_" + std::to_string(instance.d_g.id(instance.d_g.source(a))) +
                "," + std::to_string(instance.d_g.id(instance.d_g.target(a))));
    }

    // Declare flow variables.
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        for (ArcIt a(instance.d_g); a != INVALID; ++a) {
            (*f[scenarioId])[a] = model.addVar(
                0.0, instance.capacity, 0.0, GRB_CONTINUOUS,
                "f^" + std::to_string(scenarioId) + "_" +
                    std::to_string(instance.d_g.id(instance.d_g.source(a))) +
                    "," +
                    std::to_string(instance.d_g.id(instance.d_g.target(a))));
        }

        for (DNodeIt v(instance.d_g); v != INVALID; ++v) {
            if (instance.d_g.id(v) == instance.depot) {
                continue;
            }

            (*z[scenarioId])[v] =
                model.addVar(0.0, instance.capacity, 0.0, GRB_CONTINUOUS,
                             "z^" + std::to_string(scenarioId) + "_" +
                                 std::to_string(instance.d_g.id(v)));
        }
    }

    // Declare recourse variables.
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        for (DNodeIt v(instance.d_g); v != INVALID; ++v) {
            if (instance.d_g.id(v) == instance.depot) {
                continue;
            }

            (*y[scenarioId])[v] =
                model.addVar(0.0, 10.0,
                             (1.0 / float(instance.nScenarios)) *
                                 instance.getArcRecourseCost(v),
                             GRB_CONTINUOUS,
                             "y^" + std::to_string(scenarioId) + "_" +
                                 std::to_string(instance.d_g.id(v)));
        }
    }

    // Flow upper bound constraints.
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        for (ArcIt a(instance.d_g); a != INVALID; ++a) {
            model.addConstr((*f[scenarioId])[a] <= instance.capacity * x[a]);
        }

        for (DNodeIt v(instance.d_g); v != INVALID; ++v) {
            if (instance.d_g.id(v) == instance.depot) {
                continue;
            }

            model.addConstr((*z[scenarioId])[v] <=
                            instance.capacity * (*y[scenarioId])[v]);
        }
    }

    // Incidence constraints.
    for (DNodeIt v(instance.d_g); v != INVALID; ++v) {
        double RHS = (instance.d_g.id(v) == instance.depot) ? instance.k : 1.0;
        GRBLinExpr exprIn = 0.0;
        GRBLinExpr exprOut = 0.0;

        for (InArcIt a(instance.d_g, v); a != INVALID; ++a) {
            exprIn += x[a];
        }

        for (OutArcIt a(instance.d_g, v); a != INVALID; ++a) {
            exprOut += x[a];
        }

        model.addConstr(exprIn == RHS);
        model.addConstr(exprOut == RHS);
    }

    // Flow conservation.
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        for (DNodeIt v(instance.d_g); v != INVALID; ++v) {
            if (instance.d_g.id(v) == instance.depot) {
                continue;
            }

            GRBLinExpr expr = 0.0;

            for (InArcIt a(instance.d_g, v); a != INVALID; ++a) {
                expr -= (*f[scenarioId])[a];
            }

            for (OutArcIt a(instance.d_g, v); a != INVALID; ++a) {
                expr += (*f[scenarioId])[a];
            }

            model.addConstr(
                (*z[scenarioId])[v] + expr ==
                instance.scenariosMatrix[instance.d_g.id(v)][scenarioId]);
        }
    }
}

bool FlowModel::solve(SVRPSolution &solution) {
    // First we solve the root LP.
    solveRootLP();
    solution.rootBound = model.get(GRB_DoubleAttr_ObjVal);

    for (ArcIt a(instance.d_g); a != INVALID; ++a) {
        if (x[a].get(GRB_DoubleAttr_X) > 0.001) {
            std::cout << x[a].get(GRB_StringAttr_VarName) << " = "
                      << x[a].get(GRB_DoubleAttr_X) << std::endl;
        }
    }

    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        for (DNodeIt v(instance.d_g); v != INVALID; ++v) {
            if (instance.d_g.id(v) != instance.depot) {
                if ((*y[scenarioId])[v].get(GRB_DoubleAttr_X) > 0.001) {
                    std::cout << (*y[scenarioId])[v].get(GRB_StringAttr_VarName)
                              << " = "
                              << (*y[scenarioId])[v].get(GRB_DoubleAttr_X)
                              << std::endl;
                }
            }
        }
    }

    // Change variables type.
    // for (ArcIt a(instance.d_g); a != INVALID; ++a) {
    //     x[a].set(GRB_CharAttr_VType, GRB_BINARY);
    // }

    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        for (DNodeIt v(instance.d_g); v != INVALID; ++v) {
            if (instance.d_g.id(v) != instance.depot) {
                (*y[scenarioId])[v].set(GRB_CharAttr_VType, GRB_INTEGER);
            }
        }
    }

    // Some other parameters.
    model.getEnv().set(GRB_DoubleParam_MIPGap, 1e-6);
    model.getEnv().set(GRB_DoubleParam_TimeLimit, params.timeLimit);
    model.getEnv().set(GRB_IntParam_Threads, 1);

    // Solve.
    model.update();
    model.optimize();

    // Construct solution.
    setSolution(solution);

    return true;
}

void FlowModel::solveRootLP() {
    ArcValueMap xValue(instance.d_g, 0.0);
    cvrpsepSeparator->setMaxNoOfCapCuts(20);

    bool separatedCut = true;
    while (separatedCut) {
        separatedCut = false;
        model.update();
        model.optimize();

        // Get current solution.
        for (ArcIt a(instance.d_g); a != INVALID; ++a) {
            xValue[a] = x[a].get(GRB_DoubleAttr_X);
        }

        // Separate cuts.
        std::vector<CutData> separatedCuts;
        cvrpsepSeparator->separateCVRPSEPCuts(xValue, separatedCuts);

        // Add cuts.
        if (addSeparatedCuts(separatedCuts)) {
            separatedCut = true;
        }
    }
}

bool FlowModel::addSeparatedCuts(std::vector<CutData> &separatedCuts) {
    bool addedCut = false;

    for (CutData &cutData : separatedCuts) {
        double sign = (cutData.sense == GRB_LESS_EQUAL) ? 1 : -1;
        if (sign * (cutData.LHS - cutData.RHS) <= 1e-6) {
            continue;
        }

        GRBLinExpr expr = 0.0;
        for (std::pair<Edge, double> &p : cutData.edgePairs) {
            Edge e = p.first;
            Arc a1 = findArc(
                instance.d_g,
                instance.d_g.nodeFromId(instance.g.id(instance.g.u(e))),
                instance.d_g.nodeFromId(instance.g.id(instance.g.v(e))));
            Arc a2 = findArc(
                instance.d_g,
                instance.d_g.nodeFromId(instance.g.id(instance.g.v(e))),
                instance.d_g.nodeFromId(instance.g.id(instance.g.u(e))));
            expr += p.second * x[a1] + p.second * x[a2];
        }

        // Add constraint.
        GRBConstr addedConstraint =
            model.addConstr(expr, cutData.sense, cutData.RHS, cutData.name);
        addedConstraint.set(GRB_IntAttr_Lazy, 3);
        addedCut = true;
    }

    return addedCut;
}

void FlowModel::setSolution(SVRPSolution &solution) {
    solution.lowerBound = model.get(GRB_DoubleAttr_ObjBound);
    solution.upperBound = model.get(GRB_DoubleAttr_ObjVal);
    solution.rootBound = std::max(solution.rootBound, callback->rootBound);
    solution.nodesExplored = model.get(GRB_DoubleAttr_NodeCount);
    solution.gap = model.get(GRB_DoubleAttr_MIPGap);
    solution.cvrpCuts = cvrpsepSeparator->nCuts;
    solution.deterministicCost = 0.0;
    solution.recourseCost = 0.0;
    solution.modelRecourseCost = 0.0;

    // Get back the solution.
    if (solution.upperBound <= 1e6) {
        ArcIntMap arcCount(instance.d_g);
        for (ArcIt a(instance.d_g); a != INVALID; ++a) {
            if (x[a].get(GRB_DoubleAttr_X) > 0.001) {
                arcCount[a] =
                    static_cast<int>(std::round(x[a].get(GRB_DoubleAttr_X)));
                std::cout << x[a].get(GRB_StringAttr_VarName) << " = "
                          << x[a].get(GRB_DoubleAttr_X) << std::endl;
                solution.deterministicCost +=
                    arcCount[a] * instance.d_weight[a];
            }
        }
        solution.setRoutesFromSolution(instance, arcCount);

        for (int scenarioId = 0; scenarioId < instance.nScenarios;
             scenarioId++) {
            for (DNodeIt v(instance.d_g); v != INVALID; ++v) {
                if (instance.d_g.id(v) != instance.depot) {
                    if ((*y[scenarioId])[v].get(GRB_DoubleAttr_X) > 0.001) {
                        std::cout
                            << (*y[scenarioId])[v].get(GRB_StringAttr_VarName)
                            << " = "
                            << (*y[scenarioId])[v].get(GRB_DoubleAttr_X)
                            << std::endl;
                    }
                    solution.modelRecourseCost +=
                        (1.0 / static_cast<double>(instance.nScenarios)) *
                        instance.getArcRecourseCost(v) *
                        (*y[scenarioId])[v].get(GRB_DoubleAttr_X);
                }
            }
        }

        solution.cost = solution.deterministicCost + solution.modelRecourseCost;
        solution.recourseCost = instance.getRecourseCost(solution.routes);

        std::cout << "Deterministic cost " << solution.deterministicCost
                  << endl;
        std::cout << "Model recourse cost " << solution.modelRecourseCost
                  << endl;
        std::cout << "Manually computed recourse cost " << solution.recourseCost
                  << endl;
        std::cout << "Root bound: " << solution.rootBound << std::endl;
    } else {
        std::cout << "Didn’t find optimal solution." << std::endl;
    }
}