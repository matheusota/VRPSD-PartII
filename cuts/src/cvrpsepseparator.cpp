#include "cvrpsepseparator.h"

CVRPSEPSeparator::CVRPSEPSeparator(const SVRPInstance &instance,
                                   const Params &params)
    : instance(instance), params(params), Demand(instance.n + 1, 0.0),
      NoOfCustomers(instance.n - 1), MaxNoOfCapCuts(10), MaxNoOfMStarCuts(30),
      MaxNoOfFCICuts(10), MaxNoOfCombCuts(20), MaxNoOfHypoCuts(10),
      CAP(instance.capacity), EpsForIntegrality(1e-6),
      MaxNoOfFCITreeNodes(20000) {
    // populate Demand vector
    for (NodeIt v(instance.g); v != INVALID; ++v) {
        int id = instance.g.id(v);
        if (id != 0) {
            Demand[id] = instance.demand[v];
            demandSum += Demand[id];
        }
    }

    QMin = demandSum - (instance.k - 1) * instance.capacity;
    if (params.basicVRPSD) {
        QMin = 1;
    }
}

void CVRPSEPSeparator::setMaxNoOfCapCuts(int maxNoOfCuts) {
    MaxNoOfCapCuts = maxNoOfCuts;
}

// CVRPSEP uses `N` to refer to the depot.
int CVRPSEPSeparator::checkForDepot(int i) {
    if (i == instance.n) {
        return 0;
    } else {
        return i;
    }
}

// Add capacity cuts of the form x(S:S) <= (|S| - k(S)).
int CVRPSEPSeparator::addCapacityCuts(const EdgeValueMap &xValue,
                                      const std::vector<int> &customers,
                                      int nVehicles,
                                      std::vector<CutData> &separatedCuts) {
    separatedCuts.emplace_back(CutData());
    CutData &cutData = separatedCuts.back();
    cutData.RHS = customers.size() - nVehicles;
    cutData.sense = '<';
    cutData.name = "CVRPSEP";
    cutData.LHS = 0;

    for (size_t j = 0; j < customers.size(); j++) {
        cutData.customers.push_back(customers[j]);

        for (size_t k = j + 1; k < customers.size(); k++) {
            Edge e = findEdge(instance.g, instance.g.nodeFromId(customers[j]),
                              instance.g.nodeFromId(customers[k]));
            assert(e != INVALID);
            cutData.edgePairs.push_back(std::make_pair(e, 1.0));
            cutData.LHS += xValue[e];
        }
    }

    int nSepCuts = 1;
    if (cutData.LHS <= cutData.RHS + EpsForIntegrality) {
        nSepCuts = 0;
        separatedCuts.pop_back();
    }

    nCuts += nSepCuts;
    cutData.violation = cutData.LHS - cutData.RHS;
    return nSepCuts;
}

void CVRPSEPSeparator::getEdgeMatrix(const EdgeValueMap &xValue, int &nedges,
                                     vector<vector<double>> &edgeMatrix) {
    for (EdgeIt e(instance.g); e != INVALID; ++e) {
        if (xValue[e] >= EpsForIntegrality) {
            int idu = instance.g.id(instance.g.u(e));
            int idv = instance.g.id(instance.g.v(e));

            if (idv < idu) {
                std::swap(idu, idv);
            }

            edgeMatrix[idu][idv] = xValue[e];
            nedges++;
        }
    }
}

// This maps from the directed to the undirected case.
int CVRPSEPSeparator::separateCVRPSEPCuts(const ArcValueMap &xValue,
                                          std::vector<CutData> &separatedCuts) {
    EdgeValueMap edgeValue(instance.g, 0.0);
    for (ArcIt a(instance.d_g); a != INVALID; ++a) {
        edgeValue[instance.getEdgeFromArc(a)] += xValue[a];
    }

    return separateCVRPSEPCuts(edgeValue, separatedCuts);
}

int CVRPSEPSeparator::separateCVRPSEPCuts(const EdgeValueMap &xValue,
                                          std::vector<CutData> &separatedCuts) {
    auto started = chrono::high_resolution_clock::now();

    // Count number of edges x_e > 0 and set their values to the map.
    int nedges = 0;
    vector<vector<double>> edgeMatrix(instance.n,
                                      std::vector<double>(instance.n, 0.0));
    getEdgeMatrix(xValue, nedges, edgeMatrix);

    // Populate EdgeTail, EdgeHead and EdgeX.
    std::vector<int> EdgeTail(nedges + 1), EdgeHead(nedges + 1);
    std::vector<double> EdgeX(nedges + 1);

    int k = 1;
    for (int i = 0; i < instance.n; i++) {
        for (int j = i + 1; j < instance.n; j++) {
            if (edgeMatrix[i][j] >= EpsForIntegrality) {
                int u = (i == 0) ? instance.n : i;
                int v = j;

                EdgeTail[k] = u;
                EdgeHead[k] = v;
                EdgeX[k] = edgeMatrix[i][j];
                k++;
            }
        }
    }

    // Call CVRPSEP to separate capacity cuts.
    CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;
    CMGR_CreateCMgr(&MyCutsCMP, 100);
    CMGR_CreateCMgr(&MyOldCutsCMP, 100);
    MaxCapViolation = 0;
    double newCap = (params.basicVRPSD) ? demandSum : CAP;
    CAPSEP_SeparateCapCuts(NoOfCustomers, Demand.data(), newCap, nedges,
                           EdgeTail.data(), EdgeHead.data(), EdgeX.data(),
                           MyOldCutsCMP, MaxNoOfCapCuts, EpsForIntegrality,
                           EpsForIntegrality, &IntegerAndFeasible,
                           &MaxCapViolation, MyCutsCMP);

    // Get S and k(S) from CVRPSEP.
    int nSepCuts = 0;
    customerSets.clear();
    nVehicles.clear();

    // No cuts found.
    if (IntegerAndFeasible || MyCutsCMP->Size == 0) {
        CMGR_FreeMemCMgr(&MyCutsCMP);
        CMGR_FreeMemCMgr(&MyOldCutsCMP);

        auto done = chrono::high_resolution_clock::now();
        time += static_cast<double>(
                    chrono::duration_cast<chrono::microseconds>(done - started)
                        .count()) /
                1e6;

        return 0;
    }

    for (int i = 0; i < MyCutsCMP->Size; i++) {
        if (MyCutsCMP->CPL[i]->CType == CMGR_CT_CAP) {
            assert(MyCutsCMP->CPL[i]->IntListSize >= 1);
            std::vector<int> customers;
            for (int j = 1; j <= MyCutsCMP->CPL[i]->IntListSize; j++) {
                int id = MyCutsCMP->CPL[i]->IntList[j];
                assert(id > 0 && id < instance.n);
                customers.push_back(checkForDepot(id));
            }

            int rhs = -MyCutsCMP->CPL[i]->RHS +
                      static_cast<int>(
                          customers.size()); // CVRPSEP returns |S| - k(S)
            rhs = std::max(rhs, 1);

            nSepCuts += addCapacityCuts(xValue, customers, rhs, separatedCuts);
            customerSets.push_back(customers);
            nVehicles.push_back(rhs);
        }
    }

    CMGR_FreeMemCMgr(&MyCutsCMP);
    CMGR_FreeMemCMgr(&MyOldCutsCMP);

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return nSepCuts;
}

bool CVRPSEPSeparator::updateCustomerSetsForScenario(const EdgeValueMap &xValue,
                                                     int scenarioId) {
    auto started = chrono::high_resolution_clock::now();

    // Set demands according to scenario.
    std::vector<double> scenarioDemands(NoOfCustomers + 2, 0);
    for (int i = 0; i < instance.n; i++) {
        scenarioDemands[i] =
            static_cast<double>(instance.scenariosMatrix[i][scenarioId]);
    }

    // Create empty constraint structures.
    CnstrMgrPointer scenarioCutsCMP, oldScenarioCutsCMP;
    CMGR_CreateCMgr(&scenarioCutsCMP, 100);
    CMGR_CreateCMgr(&oldScenarioCutsCMP, 100);

    // Count number of edges x_e > 0 and set their values to the map.
    int nedges = 0;
    vector<vector<double>> edgeMatrix(instance.n,
                                      std::vector<double>(instance.n, 0.0));
    getEdgeMatrix(xValue, nedges, edgeMatrix);

    // Populate EdgeTail, EdgeHead and EdgeX.
    std::vector<int> EdgeTail(nedges + 1), EdgeHead(nedges + 1);
    std::vector<double> EdgeX(nedges + 1);

    int k = 1;
    for (int i = 0; i < instance.n; i++) {
        for (int j = i + 1; j < instance.n; j++) {
            if (edgeMatrix[i][j] >= EpsForIntegrality) {
                int u = (i == 0) ? instance.n : i;
                int v = j;

                EdgeTail[k] = u;
                EdgeHead[k] = v;
                EdgeX[k] = edgeMatrix[i][j];
                k++;
            }
        }
    }

    // Call CVRPSEP to separate capacity cuts for the scenario.
    MaxCapViolation = 0;
    CAPSEP_SeparateCapCuts(
        NoOfCustomers, scenarioDemands.data(), CAP, nedges, EdgeTail.data(),
        EdgeHead.data(), EdgeX.data(), oldScenarioCutsCMP, MaxNoOfCapCuts,
        EpsForIntegrality, EpsForIntegrality, &IntegerAndFeasible,
        &MaxCapViolation, scenarioCutsCMP);

    // Get S and k(S) from CVRPSEP.
    customerSets.clear();
    nVehicles.clear();

    // No cuts found.
    if (IntegerAndFeasible || scenarioCutsCMP->Size == 0) {
        CMGR_FreeMemCMgr(&scenarioCutsCMP);
        CMGR_FreeMemCMgr(&oldScenarioCutsCMP);

        auto done = chrono::high_resolution_clock::now();
        time += static_cast<double>(
                    chrono::duration_cast<chrono::microseconds>(done - started)
                        .count()) /
                1e6;

        return false;
    }

    for (int i = 0; i < scenarioCutsCMP->Size; i++) {
        if (scenarioCutsCMP->CPL[i]->CType == CMGR_CT_CAP) {
            std::vector<int> customers;
            for (int j = 1; j <= scenarioCutsCMP->CPL[i]->IntListSize; j++) {
                int id = scenarioCutsCMP->CPL[i]->IntList[j];
                assert(id > 0 && id < instance.n);
                customers.push_back(id);
            }

            customerSets.push_back(customers);
            nVehicles.push_back(
                std::max(static_cast<int>(-scenarioCutsCMP->CPL[i]->RHS) +
                             static_cast<int>(customers.size()),
                         1)); // CVRPSEP returns |S| - k(S)
        }
    }

    CMGR_FreeMemCMgr(&scenarioCutsCMP);
    CMGR_FreeMemCMgr(&oldScenarioCutsCMP);

    auto done = chrono::high_resolution_clock::now();
    time += static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(done - started)
                    .count()) /
            1e6;

    return true;
}

// void CVRPSEPSeparator::checkMyOldCuts() const {
//     for (int i = 0; i < MyOldCutsCMP->Size; i++) {
//         if (MyOldCutsCMP->CPL[i]->CType == CMGR_CT_CAP) {
//             // if (i == 2198) {
//             //     std::cout << "CUT " << i << std::endl;
//             // }
//             for (int j = 1; j <= MyOldCutsCMP->CPL[i]->IntListSize; j++) {
//                 int id = MyOldCutsCMP->CPL[i]->IntList[j];
//                 assert(id >= 1 && id < instance.n);
//                 // if (i == 2198) {
//                 //     std::cout << id << std::endl;
//                 // }
//             }
//         }
//     }
// }