#include "persistence.h"

BoundaryMatrix::BoundaryMatrix() {
}

BoundaryMatrix::~BoundaryMatrix() {
}

BoundaryMatrix::BoundaryMatrix(std::vector<int> &_h, std::vector< std::list<int> > &_d) {
    this->header = _h;
    this->data = _d;
}

void BoundaryMatrix::print() {
    for(auto& col : data) {
        for(auto& row : col) {
            std::cout << row << " ";
        }
        std::cout << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////

PersistentHomology::PersistentHomology()  {
}

PersistentHomology::~PersistentHomology()  {
}

BoundaryMatrix PersistentHomology::compute_matrix(const SimplicialComplex &sc) {
    std::map<std::string, int> mapping;
    return PersistentHomology::compute_matrix(sc, mapping);
}

BoundaryMatrix PersistentHomology::compute_matrix(
    const SimplicialComplex &sc,
    std::map<std::string, int> &simplex_mapping
){
	int filtration_size = sc.allSimplicis.size();

	// construct mapping between simplices and their IDs

    // if simplex_mapping is empty, it means this is not a subcomplex. Then we
    // just build the mapping for the whole complex.
    if(simplex_mapping.size() == 0) {
        for(int i = 0; i < filtration_size; i++)
            simplex_mapping[sc.allSimplicis[i].id()] = i+1;
    }

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
		Simplex simplex = sc.allSimplicis[i];

        header.push_back(simplex_mapping[simplex.id()]);

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

    BoundaryMatrix bm(header, reduction);

    // reduce the boundary matrix
    reduce_matrix2(bm);

    return bm;
}

void PersistentHomology::reduce_matrix(BoundaryMatrix &bm) {
	std::cout << "doing reduction..." << std::endl;
	ComputationTimer persistence_timer("persistence computation time");
	persistence_timer.start();

    int filtration_size = bm.data.size()-1; // minus the empty simplex
    std::vector< std::list<int> > &reduction = bm.data;

	// initialize death cycle reference - nothing there yet, so just give it all -1
    std::vector<int> death_cycle_ref(filtration_size+1);
	for(int i = 0; i < filtration_size+1; i++)
		death_cycle_ref[i] = -1;

	// perform reduction
    for(int i = 0; i < filtration_size; i++)  {
        int idx = i+1;
        std::cout << "idx: " << idx << std::endl;

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
}


std::list<int> PersistentHomology::reduce_column(std::list<int> &left_list, std::list<int> &right_list) {
    std::list<int> result;

    std::vector<int> left{std::make_move_iterator(std::begin(left_list)), std::make_move_iterator(std::end(left_list))};
    std::vector<int> right{std::make_move_iterator(std::begin(right_list)), std::make_move_iterator(std::end(right_list))};

    int i = 0;
    int j = 0;
    int ll = left.size();
    int lr = right.size();
    while(i < ll && j < lr) {
        if(left[i] < right[j]) {
            result.push_back(left[i]);
            i += 1;
        } else if(left[i] > right[j]) {
            result.push_back(right[j]);
            j += 1;
        } else {
            i += 1;
            j += 1;
        }
    }
    return result;
}

void PersistentHomology::reduce_matrix2(BoundaryMatrix &bm) {
	std::cout << "doing reduction2..." << std::endl;
	ComputationTimer persistence_timer("persistence computation time");
	persistence_timer.start();

    int filtration_size = bm.data.size()-1; // minus the empty simplex
    std::vector< std::list<int> > &m = bm.data;
    std::map<int, int> low_to_col; // from low[i] to col indices

    int i = -1;
    for(auto &col : m) {
        i ++;
        if(col.size() == 0) continue;

        int mx = col.back();
        while(low_to_col.find(mx) != low_to_col.end()) {
            auto& col_to_reduce = m[low_to_col[mx]];
            //TODO More efficient: we can modify col in place
            auto new_col = reduce_column(col_to_reduce, col);
            col = new_col;
            if(col.size() == 0) {
                mx = -1;
                break;
            }
            mx = col.back();
        }

        if(mx == -1) {
            //if(low_to_col.find(bm.header[i]) != low_to_col.end()) {
                //low_to_col.erase();
            //}
        } else {
            low_to_col[mx] = i;
        }
    }

	persistence_timer.end();
	persistence_timer.dump_time();
}

BoundaryMatrix PersistentHomology::compute_matrix( Cover &cover ) {
    std::cout << "calculating reduction matrix for cover\n";

    std::vector<BoundaryMatrix> rm_vec; // reduced matrices vector
    for(int i = 0; i < cover.subComplexSize(); ++ i) {
        std::cout << "subcomplex: " << cover.IDs[i] << std::endl;
        //cover.subComplexes[cover.IDs[i]].print();

        BoundaryMatrix bm = PersistentHomology::compute_matrix(
            cover.subComplexes[cover.IDs[i]],
            cover.SimplexIDMap
        );
        rm_vec.push_back(bm);
    }

    std::cout << "gluing...\n";
    // TODO glue them together
    // Assumption: each subcomplex is already sorted
    // 1. reorder each column
    // 2. reduce the new matrix
    // Since we only have a small number of subcomplexes, just use linear scan
    // to merge them. Optimally, we can use a heap.
    std::vector<int> merged_header;
    std::vector< std::list<int> > merged_data;

    // FIXME Only a place holder
    for(int i = 0; i < cover.subComplexSize(); i ++) {
        for(int j = 0; j < rm_vec[i].size(); j ++) {
            merged_header.push_back(rm_vec[i].header[j]);
            merged_data.push_back(rm_vec[i].data[j]);
        }
    }

    /*
    // pointers to the head of each subcomplex vector
    // start from 1 to ignore the empty simplex for now
    std::vector<int> p(cover.subComplexSize(), 1);

    merged_header.push_back(0);
    merged_data.push_back(std::list<int>());

    while(true) {
        int minHead = 0; // which subcomplex has the min head
        for(int i = 1; i < cover.subComplexSize(); i ++) {
            if(p[i] < rm_vec[i].size() &&
               rm_vec[i].header[p[i]] < rm_vec[minHead].header[p[minHead]]) {
                minHead = i;
            }
        }
        merged_header.push_back(p[minHead]);
        merged_data.push_back(rm_vec[minHead][p[minHead]]);
        p[minHead] ++;
        if(p[minHead] == rm_vec[minHead].size()) {
            p[minHead] = -1;
        }

        int sum = 0;
        for(auto& n : p) {
            sum += n;
        }

        if(sum == -1*p.size()) {
            break;
        }
    }
    */

    std::cout << "reducing glued matrix\n";
    BoundaryMatrix bm(merged_header, merged_data);
    reduce_matrix2(bm);
    return bm;
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
