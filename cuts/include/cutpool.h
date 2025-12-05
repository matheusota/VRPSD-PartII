#ifndef CUTPOOL_H
#define CUTPOOL_H

#include "cutdata.h"
#include "mygraphlib.h"
#include <vector>

class CutPool {
  public:
    double rootBound = 0.0;
    CutPool() {};

    int addCut(CutData &cutData) {
        cuts.push_back(cutData);
        return static_cast<int>(cuts.size()) - 1;
    }

    CutData &getCut(int index) {
        assert(index >= 0 && static_cast<size_t>(index) < cuts.size());
        return cuts[index];
    }

    const std::vector<CutData> &getAllCuts() { return cuts; }

    void eraseCuts() {
        rootBound = 0.0;
        cuts.clear();
    }

  private:
    std::vector<CutData> cuts;
};
#endif // CUTPOOL_H
