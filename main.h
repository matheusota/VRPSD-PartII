#ifndef MAIN_H
#define MAIN_H

#include "filehandler.h"
#include "mygraphlib.h"
#include "myutils.h"
#include "params.h"
#include "stringhelper.h"
#include "svrpinstance.h"
#include "svrpsolution.h"
#include "graphviewer.h"
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <unordered_map>
#include <vector>
#include <chrono>
#include "noderecoursemodel.h"
#include "nodescenariosmodel.h"
#include "cutpool.h"
#include "cxxopts.h"

Params getParams(int argc, const char *argv[]);
void printSolutionDetails(const SVRPSolution &solution);

#endif
