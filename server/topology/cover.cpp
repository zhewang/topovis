#include <algorithm>

#include "cover.h"
#include "string"

Cover::Cover() {
}

Cover::Cover(const SimplicialComplex &sc, std::map<int,int> &vertex_map) {
    // Get the simplex ID mapping
    this->SimplexIDMap = sc.get_simplex_map();

    // TODO we can only keep the indices for each subcomplex to reduce memory
    // usage

    // Calculate sub-complexes
    std::map<int, std::vector<Simplex> > subcomplex_map;

    for(auto it = sc.allSimplicis.begin(); it != sc.allSimplicis.end(); it ++) {
        int dest = vertex_map[it->min_vertex()];

        if(subcomplex_map.find(dest) == subcomplex_map.end()) {
            subcomplex_map[dest] = std::vector<Simplex>();
        }
        subcomplex_map[dest].push_back(*it);
    }

    std::vector< std::vector<int> > subcomplex_IDs;
    for(auto it = subcomplex_map.begin(); it != subcomplex_map.end(); it ++) {
        subComplexes.push_back(SimplicialComplex(it->second, true));
        subcomplex_IDs.push_back(std::vector<int>());
        for(auto& e : subComplexes.back().allSimplicis) {
            subcomplex_IDs.back().push_back(SimplexIDMap[e.id()]);
        }
    }

    // calculate blowup complex
    if(subcomplex_IDs.size() > 1) {
        this->intersection = subcomplex_IDs[0];

        //sc.print();

        for(int i = 0; i < subComplexes.size(); i ++) {
            //subComplexes[IDs[i]].print();
            //for(auto & e: subcomplex_IDs[IDs[i]]) {
                //std::cout << e << " ";
            //}
            //std::cout << std::endl;

            std::vector<int> temp;
            std::set_intersection(subcomplex_IDs[i].begin(),
                                  subcomplex_IDs[i].end(),
                                  intersection.begin(),
                                  intersection.end(),
                                  std::back_inserter(temp));
            this->intersection = temp;
        }

        std::cout << "blowup complex intersection: " << std::endl;
        std::vector<Simplex> intersection_simplicis;
        int count = 0;
        for(auto& id: intersection) {
            if(count != 0) {
                std::cout << " ";
            }
            count ++;
            std::cout << "{" << sc.allSimplicis[id-1].id() << "}";

            intersection_simplicis.push_back(sc.allSimplicis[id-1]);

        }
        std::cout << std::endl;
    }

}
