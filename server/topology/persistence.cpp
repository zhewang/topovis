#include "persistence.h"
#include <algorithm>


BMCell::BMCell() {
    first = 0;
    second = 0;
}

BMCol::BMCol() {
    header = BMCell();
    faces = std::vector<BMCell>();
}

void BMCol::print() {
    std::cout << "----" << std::endl;
    header.print();
    std::cout << "...." << std::endl;
    for(auto & e: faces) {
        e.print();
    }
}

BMatrix::BMatrix(std::vector<BMCol> &_cols) {
    this->cols = _cols;
}
///////////////////////////////////////////////////////////////////////////////

PersistentHomology::PersistentHomology()  {
}

PersistentHomology::~PersistentHomology()  {
}

BMatrix PersistentHomology::compute_matrix(const SimplicialComplex &sc) {
    std::map<std::string, int> mapping;
    return PersistentHomology::compute_matrix(sc, 0, mapping);
}

BMatrix PersistentHomology::compute_matrix(
    const SimplicialComplex &sc,
    int complexID,
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
    std::vector<BMCol> cols;
    cols.resize(filtration_size+1);

	// reserve 1st entry as dummy simplex, in line with reduced persistence
	cols[0] = BMCol();

	// for each simplex, take its boundary and assign those simplices to its list
	for(int i = 0; i < filtration_size; i++)  {
		int idx = i+1;
		Simplex simplex = sc.allSimplicis[i];
        int globalID = simplex_mapping[simplex.id()];

        cols[idx] = BMCol();
        cols[idx].header = BMCell(globalID, complexID);

		// if 0-simplex, then reserve face as dummy simplex
		if(simplex.dim()==0)  {
			cols[idx].faces.push_back(BMCell(0, complexID));
			continue;
		}

		std::vector<Simplex> faces = simplex.faces();
		for(int f = 0; f < faces.size(); f++)  {
			Simplex next_face = faces[f];
			int face_id = simplex_mapping[next_face.id()];
			cols[idx].faces.push_back(BMCell(face_id, complexID));
		}
		// sort list, so we can efficiently add cycles and inspect death cycles
        std::sort(cols[idx].faces.begin(), cols[idx].faces.end());
	}

    BMatrix bm(cols);

    // reduce the boundary matrix
    reduce_matrix2(bm);

    return bm;
}

//void PersistentHomology::reduce_matrix(BoundaryMatrix &bm) {
	//std::cout << "doing reduction..." << std::endl;
	//ComputationTimer persistence_timer("persistence computation time");
	//persistence_timer.start();

    //int filtration_size = bm.data.size()-1; // minus the empty simplex
    //std::vector< std::list<int> > &reduction = bm.data;

	//// initialize death cycle reference - nothing there yet, so just give it all -1
    //std::vector<int> death_cycle_ref(filtration_size+1);
	//for(int i = 0; i < filtration_size+1; i++)
		//death_cycle_ref[i] = -1;

	//// perform reduction
    //for(int i = 0; i < filtration_size; i++)  {
        //int idx = i+1;
        //std::cout << "idx: " << idx << std::endl;

        //// until we are either definitively a birth cycle or a death cycle ...
        //int low_i = reduction[idx].back();
        //int num_chains_added = 0;

        //while(reduction[idx].size() > 0 && death_cycle_ref[low_i] != -1)  {
            //num_chains_added++;
            //// add the prior death cycle to us
            //int death_cycle_ind = death_cycle_ref[low_i];

            //PHCycle::iterator our_cycle_iter = reduction[idx].begin(), added_cycle_iter = reduction[death_cycle_ind].begin();
            //while(added_cycle_iter != reduction[death_cycle_ind].end())  {
                //if(our_cycle_iter == reduction[idx].end())  {
                    //reduction[idx].push_back(*added_cycle_iter);
                    //++added_cycle_iter;
                    //continue;
                //}
                //int sigma_1 = *our_cycle_iter, sigma_2 = *added_cycle_iter;
                //if(sigma_1 == sigma_2)  {
                    //our_cycle_iter = reduction[idx].erase(our_cycle_iter);
                    //++added_cycle_iter;
                //}
                //else if(sigma_1 < sigma_2)
                    //++our_cycle_iter;
                //else  {
                    //reduction[idx].insert(our_cycle_iter, sigma_2);
                    //++added_cycle_iter;
                //}
            //}
            //low_i = reduction[idx].back();
        //}

        //// if we are a death cycle then add us to the list, add as persistence pairings
        //if(reduction[idx].size() > 0)  {
            //death_cycle_ref[low_i] = idx;
            //// kill cycle at low_i, since it represents a birth
            //reduction[low_i] = PHCycle();
        //}
    //}

	//persistence_timer.end();
	//persistence_timer.dump_time();
//}


BMCol PersistentHomology::reduce_column(BMCol &left_col, BMCol &right_col) {
    BMCol result;
    result.header = left_col.header;

    auto & left = left_col.faces;
    auto & right = right_col.faces;

    int il = 0;
    int ir = 0;
    while(il < left.size() && ir < right.size()) {
        if(left[il] < right[ir]) {
            result.faces.push_back(left[il]);
            il += 1;
        } else if(right[ir] < left[il]) {
            result.faces.push_back(right[ir]);
            ir += 1;
        } else {
            il += 1;
            ir += 1;
        }
    }
    return result;
}

void PersistentHomology::reduce_matrix2(BMatrix &bm) {
	std::cout << "doing reduction2..." << std::endl;
	ComputationTimer persistence_timer("persistence computation time");
	persistence_timer.start();

    int filtration_size = bm.cols.size()-1; // minus the empty simplex
    std::map<BMCell, int> low_to_col; // from low[i] to col indices

    int i = -1;
    BMCell minusOne = BMCell(-1, -1);
    for(auto &col : bm.cols) {
        i ++;
        if(col.faces.size() == 0) continue;

        BMCell mx = col.faces.back();
        while(low_to_col.find(mx) != low_to_col.end()) {
            auto& col_to_reduce = bm.cols[low_to_col[mx]];
            auto new_col = reduce_column(col, col_to_reduce);
            col = new_col;
            if(col.faces.size() == 0) {
                mx = minusOne;
                break;
            }
            mx = col.faces.back();
        }

        if(mx == minusOne) {
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

BMatrix PersistentHomology::compute_matrix( Cover &cover ) {
    std::cout << "calculating reduction matrix for cover\n";

    std::vector<BMatrix> rm_vec; // reduced matrices vector
    for(int i = 0; i < cover.subComplexCount(); ++ i) {
        std::cout << "subcomplex: " << i << std::endl;
        //cover.subComplexes[cover.IDs[i]].print();

        BMatrix bm = PersistentHomology::compute_matrix(
            cover.subComplexes[i],
            i,
            cover.SimplexIDMap
        );
        rm_vec.push_back(bm);
    }

    std::cout << "gluing...\n";
    BMatrix bm = rm_vec[0]; // FIXME just place holder

    std::cout << "reducing glued matrix\n";
    reduce_matrix2(bm);
    return bm;
}

PersistenceDiagram* PersistentHomology::read_persistence_diagram
(BMatrix &reduction, SimplicialComplex &sc) {

	int filtration_size = sc.allSimplicis.size();

	std::vector< std::pair<int,int> > persistence_pairing;

    for(int i = 1; i < reduction.cols.size(); i++)  {

        // until we are either definitively a birth cycle or a death cycle ...
        if(reduction.cols[i].faces.size() == 0) {
            continue;
        }

        int low_i = reduction.cols[i].faces.back().first;

        // if we are a death cycle then add us to the list, add as persistence pairings
        if(reduction.cols[i].faces.size() > 0)  {
            if(low_i > 0)
                persistence_pairing.push_back(std::pair<int,int>(
                            low_i-1,
                            reduction.cols[i].header.first-1));
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
