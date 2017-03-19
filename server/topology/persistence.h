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

class PersistentHomology  {
	public:
		PersistentHomology();
		~PersistentHomology();

        //bool compute_matrix(std::vector<PHCycle> &reduction);
        static std::vector<PHCycle> compute_matrix( SimplicialComplex &sc );
        static PersistenceDiagram* compute_persistence(std::vector<PHCycle> &reduction, SimplicialComplex &sc);

        // TODO test, result from the whole complex should be the same from the cover with blowup complex
        // static std::vector<PHCycle> compute_matrix( Cover &cover );
};

#endif
