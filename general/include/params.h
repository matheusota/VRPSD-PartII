#ifndef PARAMS_H
#define PARAMS_H

#include <string>

enum POLICY { CLASSICAL, SCENARIO_OPTIMAL };

typedef struct structParams {
    POLICY policy = CLASSICAL;
    int timeLimit = 0;
    int k = -1;
    bool readKfromFile = false;
    bool showGraph = false;
    std::string inputFile = "";
    std::string outputFile = "";
    std::string tableName = "";
    bool paradaSet = false;
    bool sriCuts = false;
    bool partialRouteCuts = false;
    bool basicVRPSD = false;
    bool sriFlowSeparation = false;
    bool inOut = false;
} Params;
#endif // PARAMS_H
