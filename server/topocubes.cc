#include "topocubes.h"

TopoCubes::TopoCubes() {
}

TopoCubes::TopoCubes(std::string filePath, int _max_d) {
  ifstream csvfile;
  csvfile.open(filePath);

  std::string header;
  double x,y,c;
  std::getline(csvfile, header);
  while(csvfile >> x >> y >> c) {
    vector<double> p, a;
    p.push_back(x);
    p.push_back(y);
    this->points.push_back(Vector(p));

    this->vertex_map[points.size()-1] = c;

    a.push_back(c);
    this->attrs.push_back(a);
  }

  csvfile.close();
  this->max_d = _max_d;

  this->BuildCube();
}

void TopoCubes::BuildCube() {
    //auto vmap = this->get_quadtree_map(this->points);
    auto vmap = this->get_category_map(0);

    Filtration* filtration = new RipsFiltration(points, this->max_d);
    //Filtration* filtration = new SparseRipsFiltration(points, this->max_d, 1.0/3);
    filtration->build_filtration();
    this->global_complex = filtration->get_complex();

    global_compare::order_map = this->global_complex.get_simplex_map();
    auto globalIDMap = this->global_complex.get_simplex_map();

    for( auto it = vmap.begin(); it != vmap.end(); it ++) {
      auto s = this->getSubComplex(it->second);
      auto bm = PersistentHomology::compute_matrix(s, globalIDMap);
      // save bm to cubes
      this->cubes[it->first] = bm;
    }
}

SimplicialComplex TopoCubes::getSubComplex(std::vector<int> &ids) {
  std::vector<Simplex> subComplex;
  for(auto it = this->global_complex.allSimplicis.begin(); it != this->global_complex.allSimplicis.end(); it ++) {
    if( std::find(ids.begin(), ids.end(), it->min_vertex()) != ids.end() ) {
      subComplex.push_back(*it);
    }
  }
  return SimplicialComplex(subComplex, true);
}

std::map<int, std::vector<int> > TopoCubes::get_category_map(int offset) {
    std::map<int, std::vector<int> > vmap;
    for(int i = 0; i < this->points.size(); i ++) {
      int key = this->attrs[i][offset];
      if(vmap.find(key) == vmap.end()) vmap[key] = std::vector<int>();
      vmap[key].push_back(i);
    }
    return vmap;
}

std::vector<int> TopoCubes::parseQuery(json q) {
  std::vector<int> parsed;
  for (json::iterator it = q.begin(); it != q.end(); ++it) {
    parsed.push_back(*it);
  }
  return parsed;
}

json TopoCubes::queryCategories(json query) {
  auto q = parseQuery(query);
  if(q.size() == 0) {
    return "";
  }
  // build a cover
  Cover c(this->global_complex, q, this->cubes); // Cover read simplex ids from saved bm and calculate the intersection
  // fetch bm from cubes (need to assign subcomplex ID) and calculate bm for intersection
  BMatrix reduction = PersistentHomology::compute_matrix(c, this->cubes);

  //std::cout << "reading PD..." << std::endl;
  // read pd
  PersistenceDiagram pd =
    //PersistentHomology::read_persistence_diagram(TOPOCUBES[0], gSC);
    PersistentHomology::read_persistence_diagram(reduction, this->global_complex);

  pd.sort_pairs_by_persistence();

  std::stringstream result;
  for(unsigned i = 0; i < pd.num_pairs(); i++)  {
    PersistentPair pairing = pd.get_pair(i);
    //printf("%u %.7f %.7f\n", pairing.dim(), pairing.birth_time(), pairing.death_time());
    result << pairing.dim() << " "
      << pairing.birth_time() << " "
      << pairing.death_time() << "\n";
  }

  return result.str();
}

json TopoCubes::getOriginalPointCloud() {
  json pcloud;
  for(int i = 0; i < this->points.size(); i ++) {
    json point;

    for(int d = 0; d < this->points[i].dim(); d ++) {
      point.push_back(this->points[i][d]);
    }

    for(auto & e: this->attrs[i]) {
      point.push_back(e);
    }

    pcloud.push_back(point);
  }

  json data = {};
  data["schema"] = {"px", "py", "c"};
  data["pointcloud"] = pcloud;
  return data;
}

