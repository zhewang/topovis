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
    std::map<std::string, std::vector<Simplex> > subcomplex_map;

    for(auto it = sc.allSimplicis.begin(); it != sc.allSimplicis.end(); it ++) {
        //std::cout << "min vertex: " << it->min_vertex() << std::endl;
        std::string dest = std::to_string(vertex_map[it->min_vertex()]);
        //std::cout << "dest: " << dest << std::endl;

        if(subcomplex_map.find(dest) == subcomplex_map.end()) {
            //std::cout << "new dest: " << dest << std::endl;
            this->IDs.push_back(dest);
            subcomplex_map[dest] = std::vector<Simplex>();
        }
        subcomplex_map[dest].push_back(*it);
        //std::cout << it->id() << ": " << SimplexIDMap[it->id()] << std::endl;
    }

    std::map<std::string, std::vector<int> > subcomplex_IDs;
    for(auto it = subcomplex_map.begin(); it != subcomplex_map.end(); it ++) {
        subComplexes[it->first] = SimplicialComplex(it->second, true);
        subcomplex_IDs[it->first] = std::vector<int>();
        for(auto& e : subComplexes[it->first].allSimplicis) {
            subcomplex_IDs[it->first].push_back(SimplexIDMap[e.id()]);
        }
    }

    // calculate blowup complex
    std::vector<int> intersection = subcomplex_IDs[IDs[0]];

    //sc.print();

    for(int i = 0; i < this->IDs.size(); i ++) {
        //subComplexes[IDs[i]].print();
        //for(auto & e: subcomplex_IDs[IDs[i]]) {
            //std::cout << e << " ";
        //}
        //std::cout << std::endl;

        auto it = std::set_intersection(subcomplex_IDs[IDs[i]].begin(),
                                        subcomplex_IDs[IDs[i]].end(),
                                        intersection.begin(),
                                        intersection.end(),
                                        intersection.begin());
        intersection.resize(it-intersection.begin());
    }

    //blowupComplex = SimplicialComplex(intersection);
    std::cout << "blowup complex intersection: " << std::endl;
    //sc.print();
    int count = 0;
    for(auto& id: intersection) {
        if(count != 0) {
            std::cout << " ";
        }

        std::cout << "{" << sc.allSimplicis[id-1].id() << "}";
        count ++;
        //std::cout << e << ", ";
    }
    std::cout << std::endl;
    //blowupComplex.print();

}
