#ifndef PERSISTENCE_H
#define PERSISTENCE_H

#include "persistence_diagram.h"

#include "filtration.h"
#include "sparse_rips_filtration.h"
#include "rips_filtration.h"
#include <list>
#include <map>
#include "ComputationTimer.h"

typedef std::list<int> PHCycle;

class PersistentHomology  {
	public:
		PersistentHomology(Filtration* _filtration, bool _retainGenerators=false);
		~PersistentHomology();

        PHCycle* compute_matrix();
        PersistenceDiagram* compute_persistence_from_matrix(PHCycle* reduction);

	private:
		Filtration* filtration;
		int max_d;
};

#endif
