#ifndef COVER_H
#define COVER_H

#include "simplicial_complex.h"

#include<vector>
#include<map>


class Cover {
	public:
        Cover();
		Cover(SimplicialComplex &sc, std::map<int,int> &vertex_map);

        std::map<std::string, SimplicialComplex> subComplexes;
        std::map<std::string, SimplicialComplex> blowupComplexes;

        // we need this to preserve the order of simplices in subComplexes
        std::map<std::string, int> SimplexIDMap;
};

#endif
