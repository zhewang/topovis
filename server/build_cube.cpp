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

void loadCSV(std::string filePath, Points &points, std::map<int, int> &vertex_attr_map);
void BuildCube(Points &points);
std::map<std::string, std::vector<int> > get_quadtree_map(Points &points);
SimplicialComplex getSubComplex(SimplicialComplex &global, std::vector<int> &ids);


void loadCSV(std::string filePath, Points &points, std::map<int, int> &vertex_map) {
  ifstream csvfile;
  csvfile.open(filePath);

  std::string header;
  double x,y,c;
  std::getline(csvfile, header);
  while(csvfile >> x >> y >> c) {
    vector<double> p;
    p.push_back(x);
    p.push_back(y);
    points.push_back(Vector(p));
    vertex_map[points.size()-1] = c;
  }

  csvfile.close();
}

SimplicialComplex getSubComplex(SimplicialComplex &global, std::vector<int> &ids) {
  std::vector<Simplex> subComplex;
  for(auto it = global.allSimplicis.begin(); it != global.allSimplicis.end(); it ++) {
    if( std::find(ids.begin(), ids.end(), it->min_vertex()) != ids.end() ) {
      subComplex.push_back(*it);
    }
  }
  return SimplicialComplex(subComplex);
}

void BuildCube(Points &points, std::map<std::string, BMatrix> &naiveCubes) {
    auto vmap = get_quadtree_map(points);

    Filtration* filtration = new SparseRipsFiltration(points, 2, 1.0/3);
    filtration->build_filtration();
    SimplicialComplex sc = filtration->get_complex();

    auto globalIDMap = sc.get_simplex_map();

    for( auto it = vmap.begin(); it != vmap.end(); it ++) {
      auto s = getSubComplex(sc, it->second);
      auto bm = PersistentHomology::compute_matrix(s, globalIDMap);
      // save bm to cubes
      naiveCubes[it->first] = bm;
    }
}

// divide current 2D space into 4 bins
std::map<std::string, std::vector<int> > get_quadtree_map(Points &points) {
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
    std::map<std::string, std::vector<int> > vmap;
    for(int i = 0; i < 4; i ++) { vmap[std::to_string(i)] = std::vector<int>(); }
    int id = 0;
    for(auto & p : points) {
        if(p[0] <= xmid) { // left parts
            if(p[1] <= ymid) {// bottom parts
                //vmap[id] = 1;
                vmap["1"].push_back(id);
            } else { // top parts
                //vmap[id] = 0;
                vmap["0"].push_back(id);
            }
        } else { // right parts
            if(p[1] <= ymid) {// bottom parts
                //vmap[id] = 3;
                vmap["3"].push_back(id);
            } else { // top parts
                //vmap[id] = 2;
                vmap["2"].push_back(id);
            }
        }
        id ++;
    }
    return vmap;
}

int main(int argc, char* argv[]) {
  if(argc < 2) {
    std::cout << "Usage: " << argv[0] << " [data_file_path]" << std::endl;
    return 0;
  }

  Points points;
  std::map<int, int> vertex_attr_map;
  loadCSV(argv[1], points, vertex_attr_map);

  std::map<std::string, BMatrix> naiveCubes;
  std::cout << "pass 1" << std::endl;
  BuildCube(points, naiveCubes);
  return 0;
}
