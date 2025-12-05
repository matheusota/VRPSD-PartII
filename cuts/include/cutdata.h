#ifndef CUTDATA_H
#define CUTDATA_H

#include <vector>
#include <utility>
#include <string>
#include "mygraphlib.h"

// Cuts are stored as vectors of pairs. For nodes we store (id, coef), for
// edges we store (edge, coef). We store edges differently because lemon
// treat edges ids in a strange way.
struct CutData {
    std::vector<std::pair<Edge, double>> edgePairs;
    std::vector<std::pair<int, double>> nodePairs;
    std::vector<std::pair<int, double>> scenarioPairs;
    double RHS;
    double LHS;
    double violation;
    char sense;
    std::string name;
    std::vector<int> customers;
    double dual = 0.0;
};

#endif // CUTDATA_H