#include "filtration.h"

Filtration::Filtration(int _maxD)  {
	max_d = _maxD;
}

Filtration::~Filtration()  {
}

bool Filtration::build_filtration() {
  return false;
}

bool Filtration::getSubSet(
    std::set<int> &selectedVertices,
    Filtration *selectedSimplices,
    Filtration *unselectedSimplices ) {

    for(int i = 0; i < this->all_simplices.size(); i ++) {
        Simplex s = this->get_simplex(i);

        bool selected = true;
        for(int j = 0; j < s.dim()+1; j ++) {
            if ( selectedVertices.find(s.vertex(j)) == selectedVertices.end() ) {
                selected = false;
                break;
            }
        }

        if (selected) {
            selectedSimplices->addSimplex(s);
        } else {
            unselectedSimplices->addSimplex(s);
        }
    }

    return true;
}

void Filtration::print() {
    for(int i = 0; i < this->filtration_size(); i ++) {
        std::cout << this->get_simplex(i) << std::endl;
    }
}
