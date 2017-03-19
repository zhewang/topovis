#include "simplicial_complex.h"

//struct simplexComp
//{
    //simplexComp(std::map<std::string, int> _order) { simplex_order = _order;};

    //bool operator()(const Simplex &s1, const Simplex &s2) const
    //{
        //return (this->simplex_order[s1.id()] < this->simplex_order[s2.id()]);
    //}

    //std::map<std::string, int> simplex_order;
//};

SimplicialComplex::SimplicialComplex() {

}

SimplicialComplex::SimplicialComplex
(std::vector<Simplex> &faces, bool from_faces) {
    if(from_faces) {
        std::set<Simplex> face_set(faces.begin(), faces.end());

        for(auto it = faces.begin(); it != faces.end(); it ++) {
            if(it->dim() > 0) {
                std::vector<Simplex> face_to_add = it->faces();
                for(auto it_inner = face_to_add.begin();
                    it_inner != face_to_add.end(); it_inner ++) {
                    face_set.insert(*it_inner);
                }
            }
        }

        allSimplicis = std::vector<Simplex>(face_set.begin(), face_set.end());
        std::sort(allSimplicis.begin(), allSimplicis.end());
    } else {
        allSimplicis = faces;
    }
}

std::map<std::string, int> SimplicialComplex::get_simplex_map() {
    std::map<std::string, int> map;
    for(auto it = allSimplicis.begin(); it != allSimplicis.end(); it ++) {
        map[it->id()] = map.size()+1;
    }
    return map;
}
