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

#include "topocubes.h"

using json = nlohmann::json;

static const char *s_http_port = "8800";
static struct mg_serve_http_opts s_http_server_opts;

TopoCubes CUBES;

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


    //std::cout << "reading PD..." << std::endl;
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

  sendMSG(c, result.dump());
}

static void handle_query_call2(struct mg_connection *c, struct http_message *hm) {

  json query = json::parse(string(hm->body.p, hm->body.len));
  json result = CUBES.queryCategories(query);

  sendMSG(c, result.dump());
}

static void handle_query_pointcloud(struct mg_connection *c, struct http_message *hm) {
  json pcloud = CUBES.getOriginalPointCloud();
  sendMSG(c, pcloud.dump());
}

static void ev_handler(struct mg_connection *c, int ev, void *ev_data) {
  struct http_message *hm = (struct http_message *) ev_data;

  switch (ev) {
    case MG_EV_HTTP_REQUEST:
      if (mg_vcmp(&hm->uri, "/query") == 0) {
        handle_query_call(c, hm); /* Handle RESTful call */
      } else if (mg_vcmp(&hm->uri, "/query2") == 0) {
        handle_query_call2(c, hm);
      } else if (mg_vcmp(&hm->uri, "/pointcloud") == 0) {
        handle_query_pointcloud(c, hm);
      } else {
        sendMSG(c, "TopoCubes server is running!");
      }
      break;
    default:
      break;
  }
}

int main(int argc, char *argv[]) {
  if(argc < 3) {
    std::cout << "Usage: " << argv[0] << " [data_file_path] [pre-calculated_complex_file_path]" << std::endl;
    return 0;
  }

  // build the cubes
  CUBES = TopoCubes(argv[1], argv[2], 2);

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
