#include "persistence.h"

BoundaryMatrix::BoundaryMatrix() {
}

BoundaryMatrix::~BoundaryMatrix() {
}

BoundaryMatrix::BoundaryMatrix(std::vector<int> &_h, std::vector< std::list<int> > &_d) {
    this->header = _h;
    this->data = _d;
}

PersistentHomology::PersistentHomology()  {
}

PersistentHomology::~PersistentHomology()  {
}

BoundaryMatrix PersistentHomology::compute_matrix(const SimplicialComplex &sc ) {
	int filtration_size = sc.allSimplicis.size();

	// construct mapping between simplices and their IDs
	std::map<std::string,int> simplex_mapping;
	for(int i = 0; i < filtration_size; i++)
		simplex_mapping[sc.allSimplicis[i].id()] = i+1;

	// initialize reduction to boundaries - just a vector of lists
    std::vector<int> header;
    std::vector<PHCycle> reduction;
    reduction.clear();
    reduction.resize(filtration_size+1);

	// reserve 1st entry as dummy simplex, in line with reduced persistence
	reduction[0] = PHCycle();
    header.push_back(0);

	// for each simplex, take its boundary and assign those simplices to its list
	for(int i = 0; i < filtration_size; i++)  {
		int idx = i+1;
		reduction[idx] = PHCycle();
        header.push_back(idx);
		Simplex simplex = sc.allSimplicis[i];

		// if 0-simplex, then reserve face as dummy simplex
		if(simplex.dim()==0)  {
			reduction[idx].push_back(0);
			continue;
		}

		std::vector<Simplex> faces = simplex.faces();
		for(int f = 0; f < faces.size(); f++)  {
			Simplex next_face = faces[f];
			int face_id = simplex_mapping[next_face.id()];
			reduction[idx].push_back(face_id);
		}
		// sort list, so we can efficiently add cycles and inspect death cycles
		reduction[idx].sort();
	}

    // reduce the boundary matrix

	std::cout << "doing reduction..." << std::endl;
	ComputationTimer persistence_timer("persistence computation time");
	persistence_timer.start();

	// initialize death cycle reference - nothing there yet, so just give it all -1
    std::vector<int> death_cycle_ref(filtration_size+1);
	for(int i = 0; i < filtration_size+1; i++)
		death_cycle_ref[i] = -1;

	// perform reduction
    for(int i = 0; i < filtration_size; i++)  {
        int idx = i+1;

        // until we are either definitively a birth cycle or a death cycle ...
        int low_i = reduction[idx].back();
        int num_chains_added = 0;

        while(reduction[idx].size() > 0 && death_cycle_ref[low_i] != -1)  {
            num_chains_added++;
            // add the prior death cycle to us
            int death_cycle_ind = death_cycle_ref[low_i];

            PHCycle::iterator our_cycle_iter = reduction[idx].begin(), added_cycle_iter = reduction[death_cycle_ind].begin();
            while(added_cycle_iter != reduction[death_cycle_ind].end())  {
                if(our_cycle_iter == reduction[idx].end())  {
                    reduction[idx].push_back(*added_cycle_iter);
                    ++added_cycle_iter;
                    continue;
                }
                int sigma_1 = *our_cycle_iter, sigma_2 = *added_cycle_iter;
                if(sigma_1 == sigma_2)  {
                    our_cycle_iter = reduction[idx].erase(our_cycle_iter);
                    ++added_cycle_iter;
                }
                else if(sigma_1 < sigma_2)
                    ++our_cycle_iter;
                else  {
                    reduction[idx].insert(our_cycle_iter, sigma_2);
                    ++added_cycle_iter;
                }
            }
            low_i = reduction[idx].back();
        }

        // if we are a death cycle then add us to the list, add as persistence pairings
        if(reduction[idx].size() > 0)  {
            death_cycle_ref[low_i] = idx;
            // kill cycle at low_i, since it represents a birth
            reduction[low_i] = PHCycle();
        }
    }

	persistence_timer.end();
	persistence_timer.dump_time();

    return BoundaryMatrix(header, reduction);
}

BoundaryMatrix PersistentHomology::compute_matrix( Cover &cover ) {
    std::vector<BoundaryMatrix> rm_vec; // reduced matrices vector
    for(int i = 0; i < cover.subComplexSize(); ++ i) {
        BoundaryMatrix bm = PersistentHomology::compute_matrix(
            cover.subComplexes[cover.IDs[i]]
        );
        rm_vec.push_back(bm);
    }
    return BoundaryMatrix();
}

PersistenceDiagram* PersistentHomology::read_persistence_diagram
(BoundaryMatrix &reduction, SimplicialComplex &sc) {

	int filtration_size = sc.allSimplicis.size();

	std::vector< std::pair<int,int> > persistence_pairing;

    for(int i = 0; i < filtration_size; i++)  {
        int idx = i+1;

        // until we are either definitively a birth cycle or a death cycle ...
        int low_i = reduction[idx].back();

        // if we are a death cycle then add us to the list, add as persistence pairings
        if(reduction[idx].size() > 0)  {
            if(low_i > 0)
                persistence_pairing.push_back(std::pair<int,int>(low_i-1,idx-1));
        }
    }

	std::vector<PersistentPair> persistent_pairs;
	for(int i = 0; i < persistence_pairing.size(); i++)  {
		std::pair<int,int> pairing = persistence_pairing[i];
		Simplex birth_simplex = sc.allSimplicis[pairing.first], death_simplex = sc.allSimplicis[pairing.second];
		if(death_simplex.get_simplex_distance() == birth_simplex.get_simplex_distance())
			continue;

		PersistentPair persistent_pair(birth_simplex.dim(), birth_simplex.get_simplex_distance(), death_simplex.get_simplex_distance());

		persistent_pairs.push_back(persistent_pair);
	}

	PersistenceDiagram *pd = new PersistenceDiagram(persistent_pairs);
	return pd;
}
