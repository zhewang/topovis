#include <iostream>
#include <vector>
#include <map>

#include "csv.h"
#include "geometry/covertree.h"
#include "topology/persistence.h"
#include "topology/sparse_rips_filtration.h"
#include "topology/rips_filtration.h"
#include "topology/simplicial_complex.h"
#include "topology/generalized_metric_space.h"

using namespace std;

typedef std::vector<std::vector<double> > ATTRS;

void read_mesh2D(string filePath, Points &points) {
  io::CSVReader<3, io::trim_chars<>, io::no_quote_escape<' '> > in(filePath);
  in.read_header(io::ignore_extra_column, "px", "py", "c");
  double px, py, c;
  while(in.read_row(px, py, c)){
    vector<double> p, a;
    p.push_back(px);
    p.push_back(py);
    points.push_back(Vector(p));

    //a.push_back(px); // use px as a dimension
    //a.push_back(c);
    //attrs.push_back(a);
  }
}

void read_meshSDSS(string filePath, Points &points) {
  io::CSVReader<6> in(filePath);
  in.read_header(io::ignore_extra_column, "ra", "dec", "z_photoz", "x", "y", "z");
  double ra, dec, redshift, x, y, z;
  while(in.read_row(ra, dec, redshift, x, y ,z)){
    vector<double> p, a;
    p.push_back(x);
    p.push_back(y);
    p.push_back(z);
    points.push_back(Vector(p));

    //// TODO calculate c
    //double c = 1;
    //if(redshift>= 0.031 && redshift < 0.032) {
    //c = 2;
    //} else if(redshift >= 0.032 && redshift < 0.033){
    //c = 3;
    //} else {
    //c = 4;
    //}
    //vertex_map[points.size()-1] = c;

    //a.push_back(c);
    //attrs.push_back(a);
  }
}

void process(std::string filePath, int _max_d) {

  int max_d;
  Points points;
  int originalPointsSize;
  ATTRS attrs;
  std::map<int, int> vertex_map; // vertex id -> category id
  std::map<int, BMatrix> cubes;
  SimplicialComplex global_complex;
  std::map<string,int> simplex_map;
  BMatrix global_bm;

  //read_mesh2D(filePath, points);
  read_meshSDSS(filePath, points);

  max_d = _max_d;
  originalPointsSize = points.size();

  Filtration* filtration = new RipsFiltration(points, max_d);
  //Filtration* filtration = new SparseRipsFiltration(points, max_d, 1.0/3);
  filtration->build_filtration();
  global_complex = filtration->get_complex();

  for(auto & e : global_complex.allSimplicis) {
    cout << e << " " << e.get_simplex_distance() << endl;
  }
}

int main(int argc, char* argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " PointCloudFile" << endl;
    return 0;
  }
  process(argv[1], 2);
  return 0;
}
