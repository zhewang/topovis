#ifndef PERSISTENCE_H
#define PERSISTENCE_H

#include "persistence_diagram.h"

#include "filtration.h"
#include "sparse_rips_filtration.h"
#include "rips_filtration.h"
#include "ComputationTimer.h"

#include <list>
#include <map>
#include <vector>

typedef std::list<int> PHCycle;

class PersistentHomology  {
	public:
		PersistentHomology(Filtration* _filtration, bool _retainGenerators=false);
		~PersistentHomology();

        bool compute_matrix(std::vector<PHCycle> &reduction);
        PersistenceDiagram* compute_persistence(std::vector<PHCycle> &reduction);

	private:
		Filtration* filtration;
		int max_d;
};

#endif
