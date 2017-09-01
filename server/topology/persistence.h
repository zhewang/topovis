#ifndef PERSISTENCE_H
#define PERSISTENCE_H

#include "persistence_diagram.h"

#include "filtration.h"
#include "sparse_rips_filtration.h"
#include "rips_filtration.h"
#include "cover.h"
#include "simplicial_complex.h"
#include "boundary_matrix.h"

#include "ComputationTimer.h"

#include <list>
#include <map>
#include <vector>

class PersistentHomology  {
	public:
		PersistentHomology();
		~PersistentHomology();

        static BMCol reduce_column(BMCol &left, BMCol &right);
        static void reduce_matrix(BMatrix &bm);

        static BMatrix compute_matrix(const SimplicialComplex &sc, std::map<std::string, int> &sIDMap);
        static BMatrix compute_matrix( Cover &cover );
        static BMatrix compute_matrix( Cover &cover, std::map<int, BMatrix> &topocubes );

        static PersistenceDiagram read_persistence_diagram(
            BMatrix &reduction,
            SimplicialComplex &sc
        );

};

#endif
