#include "persistence.h"
#include <algorithm>


BMCell::BMCell() {
    first = 0;
    second = 0;
}
///////////////////////////////////////////////////////////////////////////////

BMCol::BMCol() {
    header = BMCell(0, 0);
    faces = std::vector<BMCell>();
}

///////////////////////////////////////////////////////////////////////////////

void BMCol::print() {
    header.print();
    std::cout << "|";
    for(auto & e: faces) {
        e.print();
        std::cout << " ";
    }
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////

BMatrix::BMatrix(std::vector<BMCol> &_cols) {
    this->cols = _cols;
}

void BMatrix::append(const BMatrix &other) {
    // TODO use move_iterator for better performance
    BMCell zero(0,0);
    for(auto& e : other.cols) {
        if(e.header == zero) continue;
        this->cols.push_back(e);
    }
}

void BMatrix::sort() {
    std::sort(cols.begin(), cols.end());
}

void BMatrix::print() {
    for(auto & c: cols) {
        c.print();
    }
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
			cols[idx].faces.push_back(BMCell(0, 0));
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

    //bm.print();
    return bm;
}

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

BMatrix PersistentHomology::compute_intersection_matrix(Cover &cover) {
    std::vector<int> &ints = cover.intersection;
    int subComplexCount = cover.subComplexCount();
    std::map<std::string, int> simplex_mapping = cover.SimplexIDMap;
    SimplicialComplex &sc = cover.globalComplex;

    std::vector<BMCol> cols;
    cols.resize(ints.size());

    for(int i = 0; i < ints.size(); i ++) {
        int globalIdx = ints[i];
		Simplex simplex = sc.allSimplicis[globalIdx-1];

        cols[i] = BMCol();
        cols[i].header = BMCell(globalIdx, -1); // use -1 to represent it's intersection

        // same simplex but different subcomplex
        for(int j = 0; j < subComplexCount; j ++) {
            cols[i].faces.push_back(BMCell(globalIdx, j+1));
        }

		if(simplex.dim()==0)  {
            cols[i].faces.push_back(BMCell(0, 0));
		} else {
            std::vector<Simplex> faces = simplex.faces();
            for(int f = 0; f < faces.size(); f++)  {
                Simplex next_face = faces[f];
                int face_id = simplex_mapping[next_face.id()];
                cols[i].faces.push_back(BMCell(face_id, -1));
            }
        }
		// sort list, so we can efficiently add cycles and inspect death cycles
        std::sort(cols[i].faces.begin(), cols[i].faces.end());
    }

    BMatrix bm(cols);
    //std::cout << "intersection bm\n";
    //bm.print();
    return bm;
}

BMatrix PersistentHomology::compute_matrix( Cover &cover ) {
    std::cout << "calculating reduction matrix for cover\n";

    std::vector<BMatrix> rm_vec; // reduced matrices vector
    for(int i = 0; i < cover.subComplexCount(); ++ i) {
        std::cout << "subcomplex: " << i+1 << std::endl;
        //cover.subComplexes[cover.IDs[i]].print();

        BMatrix bm = PersistentHomology::compute_matrix(
            cover.subComplexes[i],
            i+1,
            cover.SimplexIDMap
        );
        //bm.print();
        rm_vec.push_back(bm);
    }

    //calculate boundary matrix for intersection
    rm_vec.push_back(compute_intersection_matrix(cover));

    std::cout << "gluing...\n";
    BMatrix bm = rm_vec[0];
    if(rm_vec.size() > 1) {
        for(int i = 1; i < rm_vec.size(); i ++) {
            bm.append(rm_vec[i]);
        }
    }
    bm.sort();

    std::cout << "reducing glued matrix\n";
    reduce_matrix2(bm);
    //bm.print();
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

    //for(auto &p : persistence_pairing) {
        //std::cout << "(" << p.first << ", " << p.second << ")" << std::endl;
    //}

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
