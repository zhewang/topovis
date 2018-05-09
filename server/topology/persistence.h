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

        static BMatrix compute_matrix(
            const SimplicialComplex &sc,
            std::map<string, int> &simplex_mapping,
            BMatrix &global_bm,
            std::map<int, BMatrix> &topocubes,
            std::map<int,int> &vertex_map,
            std::vector<int> q
        );

        static BMatrix compute_bm_no_reduction(const SimplicialComplex &sc, std::map<std::string, double> &dist_map);

        static PersistenceDiagram read_persistence_diagram(
            BMatrix &reduction,
            SimplicialComplex &sc
        );

};

#endif
