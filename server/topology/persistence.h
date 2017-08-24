#ifndef PERSISTENCE_H
#define PERSISTENCE_H

#include "persistence_diagram.h"

#include "filtration.h"
#include "sparse_rips_filtration.h"
#include "rips_filtration.h"
#include "cover.h"
#include "simplicial_complex.h"

#include "ComputationTimer.h"

#include <list>
#include <map>
#include <vector>


struct BMCell {
    BMCell();
    BMCell(int f, int s) { first = f; second = s; }

    int first;
    int second;

    void print() { std::cout << "(" << first << ", " << second << ")"; }

    inline bool operator==(const BMCell& other) const{
        return this->first == other.first && this->second == other.second;
    }

    inline bool operator<(const BMCell& other) const{
        if(this->first != other.first) {
            return this->first < other.first;
        } else {
            return this->second < other.second;
        }
    }
};

struct BMCol {
    BMCol();

    BMCell header;
    std::vector<BMCell> faces;

    void print();

    inline bool operator<(const BMCol& other) const{
        return this->header < other.header;
    }
};

struct BMatrix{
    BMatrix();
    BMatrix(std::vector<BMCol> &cols);

    std::vector<BMCol> cols;

    int size() { return cols.size(); }
    void append(const BMatrix &other);
    void sort();
    void print();
};

class PersistentHomology  {
	public:
		PersistentHomology();
		~PersistentHomology();

        static BMCol reduce_column(BMCol &left, BMCol &right);
        static void reduce_matrix(BMatrix &bm);

        static BMatrix compute_matrix(const SimplicialComplex &sc, std::map<std::string, int> &sIDMap);
        static BMatrix compute_matrix( Cover &cover );

        static PersistenceDiagram read_persistence_diagram(
            BMatrix &reduction,
            SimplicialComplex &sc
        );

};

#endif
