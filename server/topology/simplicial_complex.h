#ifndef SIMPLICIAL_COMPLEX_H
#define SIMPLICIAL_COMPLEX_H

#include "simplex.h"

#include<vector>
#include<map>
#include<set>

class SimplicialComplex {
    public:
        SimplicialComplex();
        SimplicialComplex(std::vector<Simplex> &faces, std::map<std::string, int> simplex_order);

        std::map<std::string, int> get_simplex_map();

        // public members
        std::vector<Simplex> allSimplicis;
};

#endif
