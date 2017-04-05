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

typedef std::list<int> PHCycle;

struct BMCell {
    BMCell(int f, int s) { first = f; second = s; }
    int first;
    int second;

    inline bool operator<(const BMCell& other) {
        return this->first < other.first;
    }
};

struct BMCol {
    BMCell header;
    std::vector<BMCell> faces;

    inline bool operator<(const BMCol& other) {
        return this->header < other.header;
    }
};

struct BMatrix{
    BMatrix();
    BMatrix(std::vector<BMCol> &cols);
    ~BMatrix();

    std::vector<BMCol> cols;

    int size() { return cols.size(); }
};

struct BoundaryMatrix {
    BoundaryMatrix();
    BoundaryMatrix(std::vector<int> &header, std::vector< std::list<int> > &data);
    ~BoundaryMatrix();

    std::vector<int> header; //the index in the global complex for each column
    std::vector< std::list<int> > data;

    int size() { return header.size(); }
    void print();

    inline std::list<int>& operator[] (const int i) {
        return this->data[i];
    }
};

class PersistentHomology  {
	public:
		PersistentHomology();
		~PersistentHomology();

        static BoundaryMatrix compute_matrix(const SimplicialComplex &sc );
        static BoundaryMatrix compute_matrix(
            const SimplicialComplex &sc,
            std::map<std::string, int> &simplex_mapping
        );
        // TODO calculate reduced matrix given a simplex mapping
        static BoundaryMatrix compute_matrix( Cover &cover );

        static std::list<int> reduce_column(std::list<int> &left, std::list<int> &right);

        static void reduce_matrix(BoundaryMatrix &bm);
        static void reduce_matrix2(BoundaryMatrix &bm);

        static PersistenceDiagram* read_persistence_diagram(BoundaryMatrix &reduction, SimplicialComplex &sc);

};

#endif
