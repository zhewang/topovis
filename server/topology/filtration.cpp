#include "filtration.h"
#include <algorithm>

Filtration::Filtration(int _maxD)  {
	max_d = _maxD;
}

Filtration::~Filtration()  {
}

bool Filtration::build_filtration() {
  return false;
}

SimplicialComplex Filtration::get_complex() {
    if(this->all_simplices.size() == 0) {
        this->build_filtration();
    }
    return SimplicialComplex(this->all_simplices);
}

bool Filtration::binarySplit(
    std::set<int> &selectedVertices,
    Filtration *selectedComplex,
    Filtration *unselectedComplex,
    Filtration *blowupComplex) {

    std::vector<Simplex> unSelected_0;
    std::vector<Simplex> unSelected_1;
    std::vector<Simplex> unSelected_2;

    for(int i = 0; i < this->filtration_size(); i ++) {
        Simplex s = this->get_simplex(i);

        bool selected = true;
        for(int j = 0; j < s.dim()+1; j ++) {
            if ( selectedVertices.find(s.vertex(j)) == selectedVertices.end() ) {
                selected = false;
                break;
            }
        }

        if (selected) {
            selectedComplex->addSimplex(s);
        } else {
            //unselectedSimplices->addSimplex(s);
            switch(s.dim()) {
                case 0: unSelected_0.push_back(s); break;
                case 1: unSelected_1.push_back(s); break;
                case 2: unSelected_2.push_back(s); break;
                default: break;
            }
        }
    }

    // unselectedSimplices is not closed by now, make it closed
    // Algorithm: for any p-simplex (p>0), if any of its vertex is not in
    // selectedVertices, add it to the complex.
    // TODO we should justify why we do this

    // missing edges from 2-simplex
    for(auto it = unSelected_2.begin(); it != unSelected_2.end(); it ++) {
        Simplex s = *it;
        std::vector<Simplex> faces = s.faces();
        for(auto face_it = faces.begin(); face_it != faces.end(); face_it ++) {
            auto find_it = std::find(unSelected_1.begin(), unSelected_1.end(), *face_it);
            if(find_it == unSelected_1.end()) {
                unSelected_1.push_back(*face_it);
                blowupComplex->addSimplex(*face_it);
            }
        }
    }


    // missing vertices from 1-simplex
    for(auto it = unSelected_1.begin(); it != unSelected_1.end(); it ++) {
        Simplex s = *it;
        std::vector<Simplex> faces = s.faces();
        for(auto face_it = faces.begin(); face_it != faces.end(); face_it ++) {
            auto find_it = std::find(unSelected_0.begin(), unSelected_0.end(), *face_it);
            if(find_it == unSelected_0.end()) {
                unSelected_0.push_back(*face_it);
                blowupComplex->addSimplex(*face_it);
            }
        }
    }

    // TODO make blowup complex closed

    // assemble unselected complex
    for(auto it = unSelected_0.begin(); it != unSelected_0.end(); it ++) {
        unselectedComplex->addSimplex(*it);
    }
    for(auto it = unSelected_1.begin(); it != unSelected_1.end(); it ++) {
        unselectedComplex->addSimplex(*it);
    }
    for(auto it = unSelected_2.begin(); it != unSelected_2.end(); it ++) {
        unselectedComplex->addSimplex(*it);
    }


    // TODO what we added are just the blowup complex??

    return true;
}

void Filtration::print() {
    std::cout << "{ ";
    for(int i = 0; i < this->filtration_size(); i ++) {
        if(i > 0) { std::cout << ", "; }
        std::cout << "{" << this->get_simplex(i) << "}";
    }
    std::cout << " }" << std::endl;
}
