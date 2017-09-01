#include "boundary_matrix.h"
#include <algorithm>


BMCell::BMCell() {
    first = 0;
    second = 0;
}
///////////////////////////////////////////////////////////////////////////////

BMCol::BMCol() {
    header = BMCell(0, 0);
    faces = std::vector<BMCell>();
}

///////////////////////////////////////////////////////////////////////////////

void BMCol::print() {
    header.print();
    std::cout << "|";
    for(auto & e: faces) {
        e.print();
        std::cout << " ";
    }
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
BMatrix::BMatrix() {
  std::vector<BMCol> cols;
  cols.push_back(BMCol());
  this->cols = cols;
}

BMatrix::BMatrix(std::vector<BMCol> &_cols) {
    this->cols = _cols;
}

void BMatrix::append(const BMatrix &other) {
    // TODO use move_iterator for better performance
    BMCell zero(0,0);
    for(auto& e : other.cols) {
        if(e.header == zero) continue;
        this->cols.push_back(e);
    }
}

void BMatrix::sort() {
    std::sort(cols.begin(), cols.end());
}

void BMatrix::print() {
    for(auto & c: cols) {
        c.print();
    }
}

std::vector<int> BMatrix::getAllSimplicis() {
  std::vector<int> simplicis;
  for(auto &c : cols) {
    simplicis.push_back(c.header.first);
  }
  return simplicis;
}

void BMatrix::setAllComplexID(int id) {
  // don't change first col
  for(int i = 1; i < cols.size(); i ++) {
    cols[i].header.second = id;
    for(auto &face : cols[i].faces) {
      face.second = id;
    }
  }
}
