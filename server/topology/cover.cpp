#include <algorithm>
#include <iterator>

#include "cover.h"
#include "string"

Cover::Cover() {
}

Cover::Cover(const SimplicialComplex &sc, std::map<int,int> &vertex_map) {
    // Get the simplex ID mapping
    this->SimplexIDMap = sc.get_simplex_map();

    this->globalComplex = sc;

    // TODO we can only keep the indices for each subcomplex to reduce memory
    // usage

    // Calculate base sub-complexes
    std::map<int, std::vector<Simplex> > subcomplex_map;

    for(auto it = sc.allSimplicis.begin(); it != sc.allSimplicis.end(); it ++) {
        int dest = vertex_map[it->min_vertex()];

        if(subcomplex_map.find(dest) == subcomplex_map.end()) {
            subcomplex_map[dest] = std::vector<Simplex>();
        }
        subcomplex_map[dest].push_back(*it);
    }

    for(auto it = subcomplex_map.begin(); it != subcomplex_map.end(); it ++) {
        std::set<int> k;
        k.insert(it->first);

        subComplexes[k] = SimplicialComplex(it->second, true);
        subcomplex_IDs[k] = std::vector<int>();

        for(auto& e : subComplexes[k].allSimplicis) {
            subcomplex_IDs[k].push_back(SimplexIDMap[e.id()]);
        }
    }

    /*
    // calculate intersection
    auto prev_set = subcomplex_IDs;
    while( prev_set.size() > 0) {
        std::map<std::set<int>, std::vector<int> > next_set;

        for(auto & prev : prev_set) {
            for(auto & base : subcomplex_IDs) {
                std::set<int> intersection_id;
                std::set_intersection(prev.first.begin(), prev.first.end(),
                                      base.first.begin(), base.first.end(),
                                      std::inserter(intersection_id, intersection_id.end()));
                if(intersection_id == base.first) {
                    // this means base complex is already in intersection
                    continue;
                }

                if( next_set.count(intersection_id) != 0){
                    // we've already computed this
                    continue;
                }

                std::vector<int> intersection_complex;
                std::set_intersection(prev.second.begin(),
                                      prev.second.end(),
                                      base.second.begin(),
                                      base.second.end(),
                                      std::back_inserter(intersection_complex));
                if(intersection_complex.size() > 0) {
                    next_set[intersection_id] = intersection_complex;
                }
            }
        }
        subcomplex_IDs.insert(next_set.begin(), next_set.end());
        prev_set = next_set;
    }
    */
}
