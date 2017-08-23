#include <algorithm>
#include <iterator>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <string>

#include "topology/simplicial_complex.h"
#include "topology/persistence.h"
#include "topology/sparse_rips_filtration.h"

void loadCSV(std::string filePath, Points &points);
void BuildCube(Points &points);
std::map<int, int> get_quadtree_map(Points &points);
void SplitComplex(const SimplicialComplex &sc, std::map<int,int> &vertex_map);


void loadCSV(std::string filePath, Points &points) {
  ifstream csvfile;
  csvfile.open(filePath);

  std::string header;
  double x,y,c;
  std::getline(csvfile, header);
  while(csvfile >> x >> y >> c) {
    cout << "x: " << x << " ";
    cout << "y: " << y << " ";
    cout << "c: " << c << " \n";
  }

  csvfile.close();
}

void BuildCube(Points &points) {
    auto vmap = get_quadtree_map(points);

    Filtration* filtration = new SparseRipsFiltration(points, 2, 1.0/3);
    filtration->build_filtration();
    SimplicialComplex sc = filtration->get_complex();

    // build a cover
    global_compare::order_map = sc.get_simplex_map();
    Cover c(sc, vmap);
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

void SplitComplex(const SimplicialComplex &sc, std::map<int,int> &vertex_map) {
    // Get the simplex ID mapping
    auto SimplexIDMap = sc.get_simplex_map();

    // Calculate base sub-complexes
    std::map<int, std::vector<Simplex> > subcomplex_map;

    for(auto it = sc.allSimplicis.begin(); it != sc.allSimplicis.end(); it ++) {
        int dest = vertex_map[it->min_vertex()];

        if(subcomplex_map.find(dest) == subcomplex_map.end()) {
            subcomplex_map[dest] = std::vector<Simplex>();
        }
        subcomplex_map[dest].push_back(*it);
    }

    std::map<std::set<int>, std::vector<int> > subcomplex_IDs; // inlcuding intersection
    for(auto it = subcomplex_map.begin(); it != subcomplex_map.end(); it ++) {
        std::set<int> k;
        k.insert(it->first);

        //subComplexes[k] = SimplicialComplex(it->second, true);
        //subcomplex_IDs[k] = std::vector<int>();

        //for(auto& e : subComplexes[k].allSimplicis) {
            //subcomplex_IDs[k].push_back(SimplexIDMap[e.id()]);
        //}
    }
}

int main(int argc, char* argv[]) {
  Points points;
  loadCSV(argv[1], points);
  return 0;
}
