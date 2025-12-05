#ifndef PARTIALROUTE_H
#define PARTIALROUTE_H

#include "mygraphlib.h"

struct PartialRouteEntry {
    std::vector<int> vertices;
    double recourseCost;
    double xInside = 0.0;
    double xLeft = 0.0;
    double xRight = 0.0;
};

struct PartialRoute {
    std::vector<PartialRouteEntry> entries;
    bool isRoute = false;
};

#endif // PARTIALROUTE_H