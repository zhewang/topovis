#include "cover.h"
#include "string"

Cover::Cover(SimplicialComplex &sc, std::map<int,int> &vertex_map) {
    std::map<std::string, std::vector<Simplex> > subcomplex_map;

    for(auto it = sc.allSimplicis.begin(); it != sc.allSimplicis.end(); it ++) {
        std::string dest = std::to_string(it->min_vertex());

        if(subcomplex_map.find(dest) == subcomplex_map.end()) {
            subcomplex_map[dest] = std::vector<Simplex>();
        }
        subcomplex_map[dest].push_back(*it);
    }

    std::map<std::string, int> simplex_order = sc.get_simplex_map();

    for(auto it = subcomplex_map.begin(); it != subcomplex_map.end(); it ++) {
        subComplexes[it->first] = SimplicialComplex(it->second, true);
    }

}
