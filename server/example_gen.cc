#include <iostream>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <typeinfo>
#include <cstdlib>
#include <fstream>

#include "json.hpp"

#include "geometry/covertree.h"

#include "topology/persistence.h"
#include "topology/sparse_rips_filtration.h"
#include "topology/rips_filtration.h"
#include "topology/simplicial_complex.h"

using json = nlohmann::json;


Points pointcloud_two_hole() {
	Points points;
    int id = 0;
    for(double x = 1; x < 17; x++) {
        for(double y = 1; y < 17; y++ ) {
            if((x == 4 || x == 12) && y >= 3 && y <= 14) continue;
            vector<double> p;
            p.push_back(x);
            p.push_back(y);
            points.push_back(Vector(p));
        }
    }
    return points;
}

void print_points(Points points){
    std::cout << "px,py,c" << std::endl;
    for(auto & v : points) {
        std::cout << v[0] << "," << v[1] << "," << 1 << std::endl;
    }
}

void print_points(Points points, std::map<int, int> vmap){
    std::cout << "px,py,c" << std::endl;
    int id = 0;
    for(auto & v : points) {
        std::cout << v[0] << "," << v[1] << "," << vmap[id]+1 << std::endl;
        id ++;
    }
}

// divide current 2D space into 4 bins
std::map<int, int> get_quadtree_map(Points &points) {
    // get the range of data
    int xmin, xmax, ymin, ymax;
    xmin = ymin = 65535;
    xmax = ymax = -1;

    for(auto & p : points) {
        if(p[0] < xmin) xmin = p[0];
        if(p[0] > xmax) xmax = p[0];
        if(p[1] < ymin) ymin = p[1];
        if(p[1] > ymax) ymax = p[1];
    }

    int xmid = (xmin+xmax)/2.0;
    int ymid = (ymin+ymax)/2.0;
    // top-left bin: 00, top-right bin: 10, bottom-left bin: 01, bottom-right bin: 11
    std::map<int, int> vmap;
    int id = 0;
    for(auto & p : points) {
        if(p[0] <= xmid) { // left parts
            if(p[1] <= ymid) {// bottom parts
                vmap[id] = 1;
            } else { // top parts
                vmap[id] = 0;
            }
        } else { // right parts
            if(p[1] <= ymid) {// bottom parts
                vmap[id] = 3;
            } else { // top parts
                vmap[id] = 2;
            }
        }
        id ++;
    }
    return vmap;
}

void compute_persistence_homology(Points &points)
{
    std::map<int,int> vertex_map;

	int max_d = 2;

    //Filtration* filtration = new RipsFiltration(points, max_d);
    //PersistentHomology ph(full_filtration);
    Filtration* filtration = new SparseRipsFiltration(points, max_d, 1.0/3);

    // and persistenthomology use this complex to calculate ph
    filtration->build_filtration();
    SimplicialComplex sc = filtration->get_complex();
    //BMatrix reduction = PersistentHomology::compute_matrix(sc);

    // build a cover
    global_compare::order_map = sc.get_simplex_map();
    Cover c(sc, vertex_map);
    BMatrix reduction = PersistentHomology::compute_matrix(c);


    std::cout << "reading PD..." << std::endl;
    // read pd
	PersistenceDiagram pd =
        PersistentHomology::read_persistence_diagram(reduction, sc);

	pd.sort_pairs_by_persistence();

    std::stringstream result;
	for(unsigned i = 0; i < pd.num_pairs(); i++)  {
		PersistentPair pairing = pd.get_pair(i);
		//printf("%u %.7f %.7f\n", pairing.dim(), pairing.birth_time(), pairing.death_time());
        result << pairing.dim() << " "
            << pairing.birth_time() << " "
            << pairing.death_time() << "\n";
	}

}

int main(int argc, char *argv[]) {
    Points points = pointcloud_two_hole();
    std::map<int,int> vmap = get_quadtree_map(points);
    print_points(points, vmap);
    return 0;
}
