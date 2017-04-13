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

void Filtration::print() {
    std::cout << "{ ";
    for(int i = 0; i < this->filtration_size(); i ++) {
        if(i > 0) { std::cout << ", "; }
        std::cout << "{" << this->get_simplex(i) << "}";
    }
    std::cout << " }" << std::endl;
}
