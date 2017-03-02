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
            selectedSimplices->addSimplex(s);
        } else {
            unselectedSimplices->addSimplex(s);
        }
    }

    // unselectedSimplices is not closed by now, make it closed
    // Algorithm: for any p-simplex (p>0), if any of its vertex is not in
    // selectedVertices, add it to the complex. By our setting, only 0-simplex
    // are needed to make the subcomplex closed.
    std::set<int> addedVertices;
    for(int i = 0; i < unselectedSimplices->filtration_size(); i ++) {
        Simplex s = unselectedSimplices->get_simplex(i);
        if(s.dim() == 0) {
            continue;
        }

        for(int j = 0; j < s.dim()+1; j ++) {
            if ( selectedVertices.find(s.vertex(j)) != selectedVertices.end() ) {
                addedVertices.insert(s.vertex(j));
            }
        }
    }

    MetricSpace* mspace = unselectedSimplices->get_simplex(0).get_metric_space();
    for(auto it = addedVertices.begin(); it != addedVertices.end(); it ++) {
        std::vector<int> add = {*it};
        unselectedSimplices->addSimplex(Simplex(add, mspace));
    }
    // TODO sort filtration

    return true;
}

void Filtration::print() {
    for(int i = 0; i < this->filtration_size(); i ++) {
        std::cout << this->get_simplex(i) << std::endl;
    }
}
