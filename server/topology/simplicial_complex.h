#ifndef SIMPLICIAL_COMPLEX_H
#define SIMPLICIAL_COMPLEX_H

#include "simplex.h"

#include<vector>
#include<map>
#include<set>

class SimplicialComplex {
    public:
        SimplicialComplex();
        SimplicialComplex(std::vector<Simplex> &faces,
                          bool from_faces = false,
                          double** distances = NULL);
        //SimplicialComplex(std::vector<Simplex> &faces, std::map<std::string, int> simplex_order);

        std::map<std::string, int> get_simplex_map() const;
        void print() const;
        std::vector<Simplex> cofacesOf(Simplex s);
        void recalculate_distances(double** distances);

        // XXX keep in mind that after deletion the complex is not a complex anymore
        void deleteSimplex(Simplex &s);

        // public members
        std::vector<Simplex> allSimplicis;
};

#endif
