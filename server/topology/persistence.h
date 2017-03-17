#ifndef PERSISTENCE_H
#define PERSISTENCE_H

#include "persistence_diagram.h"

#include "filtration.h"
#include "sparse_rips_filtration.h"
#include "rips_filtration.h"
#include "cover.h"

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
        static std::vector<PHCycle> compute_matrix( std::vector<Simplex> &sc );
        static PersistenceDiagram* compute_persistence(std::vector<PHCycle> &reduction, std::vector<Simplex> &sc);

        // TODO
        // static std::vector<PHCycle> compute_matrix( Cover &cover );
};

#endif
