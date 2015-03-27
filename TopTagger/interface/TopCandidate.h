#ifndef TOPCANDIDATE_H
#define TOPCANDIDATE_H

#include <vector>

#include "TLorentzVector.h"

enum Flavor
{
    MERGED_TOP, MERGED_W, TRIPLET_TOP, DOUBLET_TOP, DOUBLET_W
};

class TopCandidate
{
public:
    TopCandidate(const vector<int> multiplet) : multiplet_(multiplet), rsys_(nullptr), doublets_(nullptr), candVec_(nullptr) {}
    
    const vector<int>& multiplet() { return multiplet_; }

    
private:
    std::vector<int> multiplet_;
    std::vector<int> *rsys_;
    std::vector<std::vector<int>>* doublets_;
    
    TLorentzVector* candVec_;
    Flavor flavor;
};

#endif
