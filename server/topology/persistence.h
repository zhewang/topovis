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
// TODO we need a new reduced matrix struct which has header info

struct BoundaryMatrix {
    BoundaryMatrix();
    BoundaryMatrix(std::vector<int> &header, std::vector< std::list<int> > &data);
    ~BoundaryMatrix();

    std::vector<int> header;
    std::vector< std::list<int> > data;

    inline std::list<int>& operator[] (const int i) {
        return this->data[i];
    }
};

class PersistentHomology  {
	public:
		PersistentHomology();
		~PersistentHomology();

        static BoundaryMatrix compute_matrix(const SimplicialComplex &sc );
        // TODO calculate reduced matrix given a simplex mapping
        //static std::vector<PHCycle> compute_matrix( SimplicialComplex &sc, std::map<std::string, int> simplex_mapping );
        static BoundaryMatrix compute_matrix( Cover &cover );

        static PersistenceDiagram* read_persistence_diagram(BoundaryMatrix &reduction, SimplicialComplex &sc);

};

#endif
