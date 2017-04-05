#include <algorithm>

#include "cover.h"
#include "string"

Cover::Cover() {
}

Cover::Cover(SimplicialComplex &sc, std::map<int,int> &vertex_map) {
    // Get the simplex ID mapping
    this->SimplexIDMap = sc.get_simplex_map();

    // TODO we can only keep the indices for each subcomplex to reduce memory
    // usage

    // Calculate sub-complexes
    std::map<std::string, std::vector<Simplex> > subcomplex_map;

    int scCount = 0;
    for(auto it = sc.allSimplicis.begin(); it != sc.allSimplicis.end(); it ++) {
        scCount ++;
        //std::cout << "min vertex: " << it->min_vertex() << std::endl;
        std::string dest = std::to_string(vertex_map[it->min_vertex()]);
        //std::cout << "dest: " << dest << std::endl;

        if(subcomplex_map.find(dest) == subcomplex_map.end()) {
            std::cout << "new dest: " << dest << std::endl;
            this->IDs.push_back(dest);
            subcomplex_map[dest] = std::vector<Simplex>();
        }
        subcomplex_map[dest].push_back(*it);
    }

    std::cout << "simplex count: " << scCount << std::endl;

    std::map<std::string, int> simplex_order = sc.get_simplex_map();

    for(auto it = subcomplex_map.begin(); it != subcomplex_map.end(); it ++) {
        subComplexes[it->first] = SimplicialComplex(it->second, true);
    }

    // calculate blowup complex
    std::vector<Simplex> intersection(subComplexes[IDs[0]].allSimplicis.size());

    for(int i = 0; i < this->IDs.size(); i ++) {
        subComplexes[IDs[i]].print();
        auto it = std::set_intersection(subComplexes[IDs[i]].allSimplicis.begin(),
                                        subComplexes[IDs[i]].allSimplicis.end(),
                                        intersection.begin(),
                                        intersection.end(),
                                        intersection.begin());
        intersection.resize(it-intersection.begin());
    }

    blowupComplex = SimplicialComplex(intersection);
    std::cout << "blowup complex intersection: " << std::endl;
    blowupComplex.print();

}
