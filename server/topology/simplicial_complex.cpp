#include "simplicial_complex.h"


SimplicialComplex::SimplicialComplex() {

}

SimplicialComplex::SimplicialComplex
(std::vector<Simplex> &faces, bool from_faces) {
    if(from_faces) {
        std::set<Simplex, global_compare> face_set(faces.begin(), faces.end());

        for(auto it = faces.begin(); it != faces.end(); it ++) {
            if(it->dim() > 0) {
                std::vector<Simplex> face_to_add = it->faces();
                // TODO also need to insert face's faces
                for(auto it_inner = face_to_add.begin();
                    it_inner != face_to_add.end(); it_inner ++) {
                    face_set.insert(*it_inner);
                }
            }
        }

        // TODO use move_iterator for better performance
        allSimplicis = std::vector<Simplex>(face_set.begin(), face_set.end());
        global_compare c;
        //std::sort(allSimplicis.begin(), allSimplicis.end(), c);
    } else {
        allSimplicis = faces;
    }
}

std::map<std::string, int> SimplicialComplex::get_simplex_map() const {
    std::map<std::string, int> map;
    for(auto it = allSimplicis.begin(); it != allSimplicis.end(); it ++) {
        map[it->id()] = map.size()+1;
    }
    return map;
}

void SimplicialComplex::print() const {
    std::cout << "{ ";
    for(int i = 0; i < this->allSimplicis.size(); i ++) {
        if(i > 0) { std::cout << ", "; }
        std::cout << "{" << this->allSimplicis[i] << "}";
    }
    std::cout << " }" << std::endl;
}
