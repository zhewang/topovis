#ifndef COVER_H
#define COVER_H

#include "simplicial_complex.h"

#include<vector>
#include<map>


class Cover {
	public:
        Cover();
		Cover(const SimplicialComplex &sc, std::map<int,int> &vertex_map);

        int subComplexCount() const { return this->subComplexes.size(); }

        /*************** Members ****************/
        std::vector<SimplicialComplex> subComplexes;
        SimplicialComplex intersection_complex;
        SimplicialComplex globalComplex;

        std::vector<int> intersection;

        // we need this to preserve the order of simplices in subComplexes
        std::map<std::string, int> SimplexIDMap;
};

#endif
