#include "sriseparator.h"

SRISeparator::SRISeparator(const SVRPInstance &instance, const Params &params,
                           CVRPSEPSeparator &cvrpsepSeparator)
    : instance(instance), params(params), cvrpsepSeparator(cvrpsepSeparator) {}

int SRISeparator::heuristicSeparation(const EdgeValueMap &xValue,
                                      const NodeVectorValueMap &yValue,
                                      std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    int addedCuts = 0;
    for (int i = 0; i < instance.nScenarios; i++) {
        int scenarioId = instance.sortedScenarios[i].second;

        // Call CVRPSEP for each scenario to update customer sets.
        if (cvrpsepSeparator.updateCustomerSetsForScenario(xValue,
                                                           scenarioId)) {
            const std::vector<std::vector<int>> &customerSets =
                cvrpsepSeparator.getCustomerSets();

            // Try to separate node recourse cuts for the customer sets.
            for (size_t i = 0; i < customerSets.size(); i++) {
                addedCuts += addSRIFromCustomerSet(
                    xValue, yValue, customerSets[i], separatedCuts);
            }

            if (addedCuts > 0) {
                break;
            }
        }
    }

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return addedCuts;
}

int SRISeparator::mipSeparation(
    const EdgeValueMap &xValue, const NodeVectorValueMap &yValue,
    const std::chrono::time_point<std::chrono::high_resolution_clock> &started1,
    std::vector<CutData> &separatedCuts, SRI_SEPARATION_STATUS &status) {
    auto started2 = chrono::high_resolution_clock::now();

    int addedCuts = 0;
    for (int i = 0; i < instance.nScenarios; i++) {
        int scenarioId = instance.sortedScenarios[i].second;

        std::vector<int> customers;
        SRI_SEPARATION_STATUS sepStatus =
            singleScenarioMIPModel(xValue, yValue, scenarioId, customers);
        switch (sepStatus) {
        case SRI_SEPARATION_STATUS::NO_CUTS: {
            break;
        }
        case SRI_SEPARATION_STATUS::FOUND_CUTS: {
            // Try to separate node recourse cuts for the customer sets.
            addedCuts +=
                addSRIFromCustomerSet(xValue, yValue, customers, separatedCuts);
            assert(addedCuts >= 1);
            status = SRI_SEPARATION_STATUS::FOUND_CUTS;
            break;
        }
        case SRI_SEPARATION_STATUS::FAILED: {
            if (status == SRI_SEPARATION_STATUS::NO_CUTS) {
                status = SRI_SEPARATION_STATUS::FAILED;
            }
        }
        }

        if (addedCuts >= 1) {
            break;
        }

        auto done1 = chrono::high_resolution_clock::now();
        double elapsedTime =
            static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done1 - started1)
                    .count()) /
            1e6;
        if (elapsedTime > params.timeLimit) {
            break;
        }
    }

    auto done2 = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done2 - started2)
                    .count()) /
            1e6;

    return addedCuts;
}

int SRISeparator::separateSRIFromCustomerSet(
    const EdgeValueMap &xValue, const NodeVectorValueMap &yValue,
    const std::vector<int> &customers, std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    int addedCuts =
        addSRIFromCustomerSet(xValue, yValue, customers, separatedCuts);

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return addedCuts;
}

int SRISeparator::addSRIFromCustomerSet(const EdgeValueMap &xValue,
                                        const NodeVectorValueMap &yValue,
                                        const std::vector<int> &customers,
                                        std::vector<CutData> &separatedCuts) {
    // Get x(U).
    double xInside = 0.0;
    for (size_t i = 0; i < customers.size(); i++) {
        for (size_t j = i + 1; j < customers.size(); j++) {
            Edge e = findEdge(instance.g, instance.g.nodeFromId(customers[i]),
                              instance.g.nodeFromId(customers[j]));
            xInside += xValue[e];
        }
    }

    // Get scenario bounds.
    std::vector<std::pair<int, int>> scenarioBounds;
    for (int scenarioId = 0; scenarioId < instance.nScenarios; scenarioId++) {
        double scenarioSum = 0.0;
        double customersY = 0.0;

        for (const int id : customers) {
            Node v = instance.g.nodeFromId(id);
            scenarioSum += instance.scenariosMatrix[id][scenarioId];
            customersY += yValue[v][scenarioId];
        }

        int currBound = std::ceil(scenarioSum / instance.capacity);
        // Check if y^s(U) < k^s(U) + x(U) - |U|
        if (customersY <= currBound + xInside -
                              static_cast<double>(customers.size()) -
                              EpsForIntegrality) {
            scenarioBounds.push_back({scenarioId, currBound});
        }
    }

    if (scenarioBounds.empty()) {
        return 0;
    }

    return addCut(xValue, yValue, scenarioBounds, customers, separatedCuts);
}

// Add the inequality y(U) \geq \sum_{s \in S'} (k_s(U) - x(\delta(U)) / 2).
// In fact, since 2 * x(E(U)) + x(\delta(U)) = 2 * |U|, we add the
// inequality y(U) \geq \sum_{s \in S'} (k_s(U) - |U| + x(E(U))) instead.
int SRISeparator::addCut(const EdgeValueMap &xValue,
                         const NodeVectorValueMap &yValue,
                         const std::vector<std::pair<int, int>> &scenarioBounds,
                         const std::vector<int> &customers,
                         std::vector<CutData> &separatedCuts) {
    separatedCuts.emplace_back(CutData());
    CutData &cutData = separatedCuts.back();
    cutData.RHS = 0.0;
    cutData.sense = '>';
    cutData.name = "SRI";
    cutData.LHS = 0.0;
    int nScenariosCut = static_cast<int>(scenarioBounds.size());

    for (auto &[scenarioId, bound] : scenarioBounds) {
        assert(scenarioId >= 0 && scenarioId < instance.nScenarios);
        cutData.RHS += bound - static_cast<int>(customers.size());
        cutData.scenarioPairs.push_back({scenarioId, 1.0});

        for (const int id : customers) {
            assert(id >= 1 && id < instance.n);
            Node v = instance.g.nodeFromId(id);
            cutData.LHS += yValue[v][scenarioId];
        }
    }

    for (size_t i = 0; i < customers.size(); i++) {
        for (size_t j = i + 1; j < customers.size(); j++) {
            Node u = instance.g.nodeFromId(customers[i]);
            Node v = instance.g.nodeFromId(customers[j]);
            Edge e = findEdge(instance.g, u, v);
            assert(e != INVALID);
            cutData.edgePairs.push_back({e, -nScenariosCut});
            cutData.LHS += -nScenariosCut * xValue[e];
        }
    }

    for (const int id : customers) {
        assert(id >= 1 && id < instance.n);
        cutData.nodePairs.push_back({id, 1.0});
        cutData.customers.push_back(id);
    }

    if (cutData.LHS >= cutData.RHS - EpsForIntegrality) {
        separatedCuts.pop_back();
        return 0;
    } else {
        cutData.violation = cutData.RHS - cutData.LHS;
        nCuts++;
        return 1;
    }
}

SRI_SEPARATION_STATUS SRISeparator::singleScenarioMIPModel(
    const EdgeValueMap &xValue, const NodeVectorValueMap &yValue,
    int scenarioId, std::vector<int> &customers) {
    GRBEnv env(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    GRBModel model(env);

    model.set(GRB_StringAttr_ModelName, "GurobiModel");
    model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    model.set(GRB_DoubleParam_TimeLimit, 30);
    model.set(GRB_IntParam_MIPFocus, 1);
    model.set(GRB_DoubleParam_MIPGap, 1e-3);
    model.set(GRB_IntParam_Threads, 1);
    env.set(GRB_IntParam_LogToConsole, 0);

    // Set Gurobi variables.
    GRBVar scenarioVehicles =
        model.addVar(-1.0, GRB_INFINITY, 0.0, GRB_INTEGER,
                     "scenarioVehicles_" + std::to_string(scenarioId));

    NodeGRBVarMap customersBinary(instance.g);
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) == instance.depot) {
            continue;
        }

        customersBinary[v] =
            model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                         "customersBinary_" + std::to_string(instance.g.id(v)));
    }

    EdgeGRBVarMap edgesReal(instance.g);
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        Node u = instance.g.u(e);
        Node v = instance.g.v(e);
        if (instance.g.id(u) == instance.depot ||
            instance.g.id(v) == instance.depot) {
            continue;
        }

        edgesReal[e] =
            model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS,
                         "edgesReal_" + std::to_string(instance.g.id(u)) + "," +
                             std::to_string(instance.g.id(v)));
    }
    model.update();

    // Bound constraint.
    GRBLinExpr objExpr = 1 + scenarioVehicles;

    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) == instance.depot) {
            continue;
        }

        double coef = yValue[v][scenarioId];

        objExpr -= (coef + 1) * customersBinary[v];
    }

    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        Node u = instance.g.u(e);
        Node v = instance.g.v(e);
        if (instance.g.id(u) == instance.depot ||
            instance.g.id(v) == instance.depot) {
            continue;
        }

        objExpr += xValue[e] * edgesReal[e];
    }

    model.addConstr(objExpr >= 0.1);

    // Add scenario constraints.
    GRBLinExpr loadExpr = -0.1;
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) == instance.depot) {
            continue;
        }
        double demand = instance.scenariosMatrix[instance.g.id(v)][scenarioId];
        loadExpr += demand * customersBinary[v];
    }

    model.addConstr(instance.capacity * scenarioVehicles <= loadExpr);

    // Add customer constraints.
    GRBLinExpr customersExpr = 0;
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        if (instance.g.id(v) == instance.depot) {
            continue;
        }
        customersExpr += customersBinary[v];
    }
    model.addConstr(customersExpr >= 2);

    // Add edge constraints.
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        Node u = instance.g.u(e);
        Node v = instance.g.v(e);
        if (instance.g.id(u) == instance.depot ||
            instance.g.id(v) == instance.depot) {
            continue;
        }

        model.addConstr(edgesReal[e] <= customersBinary[u]);
        model.addConstr(edgesReal[e] <= customersBinary[v]);
        model.addConstr(edgesReal[e] >=
                        customersBinary[u] + customersBinary[v] - 1.0);
    }

    model.update();
    model.optimize();

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        for (NodeIt v(instance.g); v != INVALID; ++v) {
            if (instance.g.id(v) == instance.depot) {
                continue;
            }

            if (customersBinary[v].get(GRB_DoubleAttr_X) >= 0.1) {
                customers.push_back(instance.g.id(v));
            }
        }

        return SRI_SEPARATION_STATUS::FOUND_CUTS;
    } else if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
        return SRI_SEPARATION_STATUS::NO_CUTS;
    } else {
        return SRI_SEPARATION_STATUS::FAILED;
    }
}