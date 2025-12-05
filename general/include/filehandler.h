#ifndef FILEHANDLER_H
#define FILEHANDLER_H

#include "params.h"
#include "mygraphlib.h"
#include "svrpinstance.h"
#include "svrpsolution.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <regex>
#include <algorithm>

namespace filehandler {
bool readInstance(const Params &params, SVRPInstance &instance);
bool saveSolution(const SVRPSolution &sol, std::string outputName);
bool addEntryToTable(const Params &params, const SVRPInstance &instance,
                     const SVRPSolution &solution);
} // namespace filehandler

#endif // FILEHANDLER_H
