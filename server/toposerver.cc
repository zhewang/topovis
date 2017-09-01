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

#include "mongoose.h"
#include "json.hpp"

#include "geometry/covertree.h"

#include "topology/persistence.h"
#include "topology/sparse_rips_filtration.h"
#include "topology/rips_filtration.h"
#include "topology/simplicial_complex.h"


using json = nlohmann::json;

typedef std::vector<std::vector<double> > ATTRS;

static const char *s_http_port = "8800";
static struct mg_serve_http_opts s_http_server_opts;

std::map<int, BMatrix> TOPOCUBES;
SimplicialComplex gSC;

void loadCSV(std::string filePath, Points &points, ATTRS &attrs);
void BuildCube(Points &points, ATTRS &attrs, std::map<int, BMatrix> &naiveCubes);
std::map<std::string, std::vector<int> > get_quadtree_map(Points &points);
SimplicialComplex getSubComplex(SimplicialComplex &global, std::vector<int> &ids);
std::map<int, std::vector<int> > get_category_map(Points &points, ATTRS &attrs, int offset);


void loadCSV(std::string filePath, Points &points, ATTRS &attrs) {
  ifstream csvfile;
  csvfile.open(filePath);

  std::string header;
  double x,y,c;
  std::getline(csvfile, header);
  while(csvfile >> x >> y >> c) {
    vector<double> p, a;
    p.push_back(x);
    p.push_back(y);
    points.push_back(Vector(p));

    a.push_back(c);
    attrs.push_back(a);
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
  return SimplicialComplex(subComplex, true);
}

void BuildCube(Points &points, ATTRS &attrs, std::map<int, BMatrix> &naiveCubes) {
    //auto vmap = get_quadtree_map(points);
    auto vmap = get_category_map(points, attrs, 0);

    Filtration* filtration = new SparseRipsFiltration(points, 2, 1.0/3);
    filtration->build_filtration();
    gSC = filtration->get_complex();

    global_compare::order_map = gSC.get_simplex_map();

    auto globalIDMap = gSC.get_simplex_map();

    for( auto it = vmap.begin(); it != vmap.end(); it ++) {
      auto s = getSubComplex(gSC, it->second);
      auto bm = PersistentHomology::compute_matrix(s, globalIDMap);
      // save bm to cubes
      naiveCubes[it->first] = bm;
    }
}

std::map<int, std::vector<int> > get_category_map(Points &points, ATTRS &attrs, int offset) {
    std::map<int, std::vector<int> > vmap;
    for(int i = 0; i < points.size(); i ++) {
      int key = attrs[i][offset];
      if(vmap.find(key) == vmap.end()) vmap[key] = std::vector<int>();
      vmap[key].push_back(i);
    }
    return vmap;
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

void read_points_from_json(json& data, Points& points, std::map<int, int> &vertex_map)
{
    for(auto it = data["points"].begin(); it != data["points"].end(); it ++) {
        json t = *it;
        vector<double> p;
        p.push_back(t["px"]);
        p.push_back(t["py"]);
        points.push_back(Vector(p));
        vertex_map[points.size()-1] = t["c"];
        //std::cout << points.size() << ": " << t["c"] << std::endl;
    }
}

json compute_persistence_homology(json data)
{
	Points points;
    std::map<int,int> vertex_map;
	read_points_from_json(data, points, vertex_map);
  if(points.size() == 0) {
    return "0";
  }

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

    return result.str();
}

json queryCategories(std::vector<int> q) {
  // build a cover
  Cover c(gSC, q, TOPOCUBES); // Cover read simplex ids from saved bm and calculate the intersection
  // fetch bm from cubes (need to assign subcomplex ID) and calculate bm for intersection
  BMatrix reduction = PersistentHomology::compute_matrix(c, TOPOCUBES);

  std::cout << "reading PD..." << std::endl;
  // read pd
  PersistenceDiagram pd =
    PersistentHomology::read_persistence_diagram(reduction, gSC);

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

static void sendMSG(struct mg_connection *c, std::string msg) {
  const std::string sep = "\r\n";

  std::stringstream ss;
  ss << "HTTP/1.1 200 OK"             << sep
    << "Content-Type: application/json" << sep
    << "Access-Control-Allow-Origin: *" << sep
    << "Content-Length: %d"             << sep << sep
    << "%s";

  mg_printf(c, ss.str().c_str(), (int) msg.size(), msg.c_str());
}

static void handle_query_call(struct mg_connection *c, struct http_message *hm) {

  json q = json::parse(string(hm->body.p, hm->body.len));
  json result = compute_persistence_homology(q);

  /* Send result */
  sendMSG(c, result.dump());
}

std::vector<int> parseQuery(json q) {
  std::vector<int> parsed;
  for (json::iterator it = q.begin(); it != q.end(); ++it) {
    parsed.push_back(*it);
  }
  return parsed;
}

static void handle_query_call2(struct mg_connection *c, struct http_message *hm) {

  json query = json::parse(string(hm->body.p, hm->body.len));
  std::vector<int> q = parseQuery(query);
  json result = queryCategories(q);

  /* Send result */
  sendMSG(c, result.dump());
}

static void ev_handler(struct mg_connection *c, int ev, void *ev_data) {
  struct http_message *hm = (struct http_message *) ev_data;

  switch (ev) {
    case MG_EV_HTTP_REQUEST:
      if (mg_vcmp(&hm->uri, "/query") == 0) {
        handle_query_call(c, hm); /* Handle RESTful call */
      } else if (mg_vcmp(&hm->uri, "/query2") == 0) {
        handle_query_call2(c, hm);
      } else {
        sendMSG(c, "TopoCubes server is running!");
      }
      break;
    default:
      break;
  }
}

int main(int argc, char *argv[]) {
  if(argc < 2) {
    std::cout << "Usage: " << argv[0] << " [data_file_path]" << std::endl;
    return 0;
  }

  // build the cubes
  Points points;
  ATTRS attrs;
  loadCSV(argv[1], points, attrs);
  BuildCube(points, attrs, TOPOCUBES);

  // TODO build nanocubes for "ordinary" data exploration

  // start serving
  struct mg_mgr mgr;
  struct mg_connection *c;

  mg_mgr_init(&mgr, NULL);
  c = mg_bind(&mgr, s_http_port, ev_handler);
  mg_set_protocol_http_websocket(c);

  printf("Starting server on port %s\n\n", s_http_port);

  for (;;) {
    mg_mgr_poll(&mgr, 1000);
  }


  mg_mgr_free(&mgr);

  return 0;
}
