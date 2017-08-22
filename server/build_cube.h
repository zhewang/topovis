#ifndef BUILDCUBE_H
#define BUILDCUBE_H

#include<vector>
#include<map>

#include "simplicial_complex.h"
#include "topology/persistence.h"
#include "topology/sparse_rips_filtration.h"


class BuildCube{
	public:
        BuildCube();
		BuildCube(Points &point);
        void SplitComplex(const SimplicialComplex &sc, std::map<int,int> &vertex_map);

        std::map<int, int> get_quadtree_map(Points &points);
        //int subComplexCount() const { return this->subComplexes.size(); }

        //[>************** Building Datacube ***************<]
        //std::map<int, BMatrix> buildQuadTree(
            //const SimplicialComplex &sc,
            //std::map<std::string, int> SimplexIDMap,
            //int maxLevel);
        ////void split(Const SimplicialComplex &sc, std::map<int,int> &vertex_map, int level, int maxLevel)
        ////std::map<int, SimplicialComplex> subComplexes; // only base complex
        ////std::map<int, BMatrix> subComplexBMatrix; // only base complex

        //[>************** Members ***************<]
        //SimplicialComplex globalComplex;

        //std::map<std::set<int>, std::vector<int> > subcomplex_IDs; // inlcuding intersection

        //std::vector<int> intersection;

        //// we need this to preserve the order of simplices in subComplexes
        //std::map<std::string, int> SimplexIDMap;
        //// so we can only store an integer in matrix instead of a set
        //std::map<std::set<int>, int> IntersectionIDMap;
};

#endif
