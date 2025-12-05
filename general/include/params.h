#ifndef PARAMS_H
#define PARAMS_H

#include <string>

enum ALG { BASIC, NODE, NODE_SCENARIOS };
enum POLICY { CLASSICAL, SCENARIO_OPTIMAL };

typedef struct structParams {
    ALG alg = BASIC;
    POLICY policy = CLASSICAL;
    int timeLimit = 0;
    int k = -1;
    bool readKfromFile = false;
    bool showGraph = false;
    std::string inputFile = "";
    std::string outputFile = "";
    std::string tableName = "";
    bool gendreauCuts = false;
    bool paradaSet = false;
    bool sriCuts = false;
    bool partialRouteCuts = false;
    bool mipSeparation = false;
    bool scenarioOptimalPRCuts = false;
    bool benders = false;
    bool lagrangian = false;
    bool basicVRPSD = false;
    bool sriFlowSeparation = false;
    bool inOut = false;
    bool dualSetCut = false;
} Params;
#endif // PARAMS_H
