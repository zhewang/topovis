#ifndef TOPOCUBES_H
#define TOPOCUBES_H

#include <iostream>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <typeinfo>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <boost/algorithm/string.hpp>

#include "json.hpp"

#include "geometry/covertree.h"

#include "topology/persistence.h"
#include "topology/sparse_rips_filtration.h"
#include "topology/rips_filtration.h"
#include "topology/simplicial_complex.h"
#include "topology/generalized_metric_space.h"


using json = nlohmann::json;

typedef std::vector<std::vector<double> > ATTRS;

class TopoCubes {
  public:
    TopoCubes();
    TopoCubes(std::string csvFilePath, std::string complexFilePath, int max_d);
    void BuildCube();
    void BuildCube2(); // build cube on continuous attribute

    // query API
    json getOriginalPointCloud();
    json queryCategories(json q);

  private:
    SimplicialComplex getSubComplex(std::vector<int> &ids);
    std::map<int, std::vector<int> > get_category_map(int offset /*which attibute to use*/);

    std::vector<int> vector_replace( const std::vector<int> v, int old_value, int new_value);
    void subdivision();
    bool agreeingSimplex(Simplex &s);
    int getVertexLabel(int index) { return index < originalPointsSize ? vertex_map[index] : -1; };

    std::vector<int> parseQuery(json q);

    int max_d;
    Points points;
    int originalPointsSize;
    ATTRS attrs;
    std::map<int, int> vertex_map; // vertex id -> category id
    std::map<int, BMatrix> cubes;
    SimplicialComplex global_complex;
    std::map<string,int> simplex_map;
    std::map<string, double> simplex_distances;
    BMatrix global_bm;
};

#endif
