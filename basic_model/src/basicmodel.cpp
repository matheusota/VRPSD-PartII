#include "basicmodel.h"

BasicModel::BasicModel(SVRPInstance &instance, Params &params)
    : instance(instance), params(params), model(env), x(instance.g) {

    // Set model.
    setBasicModel();

    // Setup separators.
    cvrpsepSeparator = new CVRPSEPSeparator(instance, params);
    gendreauSeparator = nullptr;

    if (params.gendreauCuts) {
        gendreauSeparator = new GendreauSeparator(instance, params);
    }

    // Set Gurobi callback.
    callback = new BasicCallback(instance, params, model, x, scenariosVars,
                                 cvrpsepSeparator, gendreauSeparator);
    model.setCallback(callback);
}

BasicModel::~BasicModel() {
    delete cvrpsepSeparator;
    delete gendreauSeparator;
    delete callback;
}

void BasicModel::setBasicModel() {
    model.set(GRB_StringAttr_ModelName, "GurobiModel");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Must set to use lazy contraints.
    model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
    model.getEnv().set(GRB_IntParam_PreCrush, 1);

    // Some other parameters.
    model.getEnv().set(GRB_DoubleParam_TimeLimit, params.timeLimit);
    model.getEnv().set(GRB_IntParam_Threads, 1);
    model.getEnv().set(GRB_DoubleParam_MIPGap, 1e-6);
    model.getEnv().set(GRB_IntParam_Cuts, 0);
    model.getEnv().set(GRB_IntParam_Presolve, 2);

    // Set x variables.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        double UB = 1.0;
        if (instance.g.id(instance.g.u(e)) == instance.depot) {
            UB = instance.upperBoundMap[instance.g.v(e)];
        } else if (instance.g.id(instance.g.v(e)) == instance.depot) {
            UB = instance.upperBoundMap[instance.g.u(e)];
        }

        x[e] = model.addVar(
            0.0, UB, instance.weight[e], GRB_INTEGER,
            "x_" + std::to_string(instance.g.id(instance.g.u(e))) + "," +
                std::to_string(instance.g.id(instance.g.v(e))));
    }

    // Set recourse variables.
    for (int i = 0; i < instance.nScenarios; i++) {
        scenariosVars.push_back(model.addVar(
            0.0, GRB_INFINITY, (1.0 / static_cast<double>(instance.nScenarios)),
            GRB_CONTINUOUS, "scen_" + std::to_string(i)));
    }

    // Incidence constraints.
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        int v_id = instance.g.id(v);
        double RHS = (v_id == instance.depot) ? 2.0 * instance.k : 2.0;
        GRBLinExpr expr = 0.0;
        for (IncEdgeIt e(instance.g, v); e != INVALID; ++e) {
            expr += x[e];
        }
        model.addConstr(expr == RHS, "incidence_" + std::to_string(v_id));
    }
    model.update();
}

bool BasicModel::solve(SVRPSolution &solution) {
    double MIPTimeLimit = params.timeLimit;
    auto started = chrono::high_resolution_clock::now();

    // Solve the root separately first.
    // Change variable types.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        x[e].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
    }
    model.update();

    // Solve root.
    // cvrpsepSeparator->setMaxNoOfCapCuts(20);
    solveRootLP(started);
    solution.rootBound = model.get(GRB_DoubleAttr_ObjVal);

    // Change variable types back to integer.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        x[e].set(GRB_CharAttr_VType, GRB_INTEGER);
    }

    auto done = chrono::high_resolution_clock::now();
    solution.rootTime =
        static_cast<double>(
            chrono::duration_cast<chrono::milliseconds>(done - started)
                .count()) /
        1000.0;
    MIPTimeLimit -= solution.rootTime;

    if (MIPTimeLimit < 0) {
        return 0;
    }

    // Solve the MIP.
    model.getEnv().set(GRB_DoubleParam_TimeLimit, MIPTimeLimit);
    // cvrpsepSeparator->setMaxNoOfCapCuts(5);
    model.update();
    model.optimize();

    // Collect solution and statistics.
    setSolution(solution);

    return 0;
}

void BasicModel::setSolution(SVRPSolution &solution) {
    solution.lowerBound = model.get(GRB_DoubleAttr_ObjBound);
    solution.upperBound = model.get(GRB_DoubleAttr_ObjVal);
    solution.nodesExplored = model.get(GRB_DoubleAttr_NodeCount);
    solution.gap = model.get(GRB_DoubleAttr_MIPGap);
    solution.rootBound = std::max(solution.rootBound, callback->rootBound);
    solution.cvrpsepCuts = cvrpsepSeparator->nCuts;
    solution.deterministicCost = 0.0;
    solution.recourseCost = 0.0;
    solution.modelRecourseCost = 0.0;

    if (gendreauSeparator != nullptr) {
        solution.gendreauCuts = gendreauSeparator->nCuts;
    }

    // Get back the solution.
    if (solution.upperBound <= 1e6) {
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            solution.solved = true;
        }

        EdgeIntMap edgeCount(instance.g);
        for (EdgeIt e(instance.g); e != INVALID; ++e) {
            Node u = instance.g.u(e);
            Node v = instance.g.v(e);

            if (x[e].get(GRB_DoubleAttr_X) > 0.001) {
                edgeCount[e] =
                    static_cast<int>(std::round(x[e].get(GRB_DoubleAttr_X)));
                std::cout << "x_"
                          << std::to_string(instance.g.id(u)) + "," +
                                 std::to_string(instance.g.id(v))
                          << " = " << x[e].get(GRB_DoubleAttr_X) << std::endl;
                solution.deterministicCost += edgeCount[e] * instance.weight[e];
            }
        }
        solution.setRoutesFromSolution(instance, edgeCount);

        for (int i = 0; i < instance.nScenarios; i++) {
            double scenarioRecourse = scenariosVars[i].get(GRB_DoubleAttr_X);
            solution.modelRecourseCost +=
                1.0 / static_cast<double>(instance.nScenarios) *
                scenarioRecourse;
        }

        solution.cost = solution.deterministicCost + solution.modelRecourseCost;
        solution.recourseCost =
            instance.classicalRecourseHelper.getRecourseCost(solution.routes);
        solution.optimalRecourseCost =
            instance.optimalRecourseHelper.getRecourseCost(solution.routes);

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

void BasicModel::solveRootLP(
    std::chrono::time_point<std::chrono::high_resolution_clock> &started) {
    double prevDualBound = 0.0;
    int noImprovCount = 0;
    bool solveLP = true;
    EdgeValueMap xValue(instance.g);
    std::vector<double> scenariosRecourse(instance.nScenarios, 0.0);
    std::vector<CutData> separatedCuts;

    // model.set(GRB_IntParam_NumericFocus, 3);
    while (solveLP && noImprovCount <= 10) {
        solveLP = false;
        model.update();
        model.optimize();

        auto done = chrono::high_resolution_clock::now();
        double elapsedTime =
            static_cast<double>(
                chrono::duration_cast<chrono::milliseconds>(done - started)
                    .count()) /
            1000.0;
        if (elapsedTime > params.timeLimit) {
            return;
        }

        // Get dual bound and update no improvement counter.
        double currDualBound = model.get(GRB_DoubleAttr_ObjBound);
        if (currDualBound - prevDualBound >= 1e-3) {
            noImprovCount = 0;
        } else {
            noImprovCount++;
        }
        prevDualBound = currDualBound;

        // Construct current solution.
        double totalRecourse = 0.0;
        getCurrentSolution(xValue, scenariosRecourse, totalRecourse);

        separatedCuts.clear();
        cvrpsepSeparator->separateCVRPSEPCuts(xValue, separatedCuts);

        if (gendreauSeparator != nullptr && separatedCuts.empty()) {
            gendreauSeparator->separateAggregatedGendreauCut(
                xValue, totalRecourse, separatedCuts);
        }

        if (addSeparatedCuts(separatedCuts)) {
            noImprovCount = 0;
            solveLP = true;
            continue;
        }
    }
}

void BasicModel::getCurrentSolution(EdgeValueMap &xValue,
                                    std::vector<double> &scenarioRecourses,
                                    double &totalRecourse) {
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        xValue[e] = x[e].get(GRB_DoubleAttr_X);
    }

    // Construct recourseValue.
    for (int s = 0; s < instance.nScenarios; s++) {
        scenarioRecourses[s] = scenariosVars[s].get(GRB_DoubleAttr_X);
        totalRecourse += scenarioRecourses[s];
    }
}

bool BasicModel::addSeparatedCuts(std::vector<CutData> &separatedCuts) {
    bool addedCut = false;

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

        // std::cout << cutData.name << std::endl;
        // std::cout << expr << " " << cutData.sense << " " << cutData.RHS
        //           << std::endl;
        GRBConstr addedConstraint =
            model.addConstr(expr, cutData.sense, cutData.RHS, cutData.name);
        addedConstraint.set(GRB_IntAttr_Lazy, 3);
        addedCut = true;
    }

    return addedCut;
}