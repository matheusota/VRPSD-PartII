#ifndef CVRPSEPSEPARATOR_H
#define CVRPSEPSEPARATOR_H
// #define SCIP_DEBUG
#include "mygraphlib.h"
#include "svrpinstance.h"
#include "params.h"
#include <cassert>
#include <lemon/list_graph.h>
#include <list>
#include <utility>
#include <ctime> // For CLOCKS_PER_SEC
#include <chrono>
#include "cnstrmgr.h"
#include "capsep.h"
#include "mstarsep.h"
#include "fcisep.h"
#include "combsep.h"
#include "htoursep.h"
#include "cnstrmgr.h"
#include "cutdata.h"

class CVRPSEPSeparator {
  public:
    int nCuts = 0;
    double time = 0.0;

    CVRPSEPSeparator(const SVRPInstance &instance, const Params &params);
    void setMaxNoOfCapCuts(int maxNoOfCuts);
    const std::vector<std::vector<int>> &getCustomerSets() const {
        return customerSets;
    };
    const std::vector<int> &getNVehicles() const { return nVehicles; };

    // These functions return the number of separated cuts.
    int separateCVRPSEPCuts(const EdgeValueMap &xValue,
                            std::vector<CutData> &separatedCuts);

    // Used for other separations that uses CVRPSEP.
    bool updateCustomerSetsForScenario(const EdgeValueMap &xValue,
                                       int scenarioId);

    // void checkMyOldCuts() const;

  private:
    const SVRPInstance &instance;
    const Params &params;

    // Parameters for the CVRPSEP package.
    double demandSum = 0;
    std::vector<double> Demand;
    int NoOfCustomers, NoOfEdges, MaxNoOfCapCuts, MaxNoOfMStarCuts,
        MaxNoOfFCICuts, MaxNoOfCombCuts, MaxNoOfHypoCuts;
    // CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;
    double CAP, EpsForIntegrality, MaxCapViolation, MaxMStarViolation,
        MaxFCIViolation, MaxCombViolation, MaxHypoViolation;
    char IntegerAndFeasible;
    int MaxNoOfFCITreeNodes;
    int QMin;

    // The separator remembers the customer sets S and the RHS' k(S) so that one
    // can use these later if desired.
    std::vector<std::vector<int>> customerSets;
    std::vector<int> nVehicles;

    int checkForDepot(int i);
    void getEdgeMatrix(const EdgeValueMap &xValue, int &nedges,
                       vector<vector<double>> &edgeMatrix);
    int addCapacityCuts(const EdgeValueMap &xValue,
                        const std::vector<int> &customers, int nVehicles,
                        std::vector<CutData> &separatedCuts);
};

#endif // CVRPSEPSEPARATOR_H
