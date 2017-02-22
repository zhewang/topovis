#include <iostream>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <typeinfo>
#include <cstdlib>

//#include <boost/random.hpp>
//#include <boost/generator_iterator.hpp>
//#include <boost/date_time/gregorian/gregorian.hpp>
//#include <boost/date_time/posix_time/posix_time.hpp>
//#include <boost/algorithm/string.hpp>

#include "mongoose.h"
#include "json.hpp"

#include "topology/persistence.h"
#include "geometry/covertree.h"
#include "topology/sparse_rips_filtration.h"
#include "topology/rips_filtration.h"

using json = nlohmann::json;


static const char *s_http_port = "8800";
static struct mg_serve_http_opts s_http_server_opts;

void read_points_from_json(json& data, Points& points)
{
    for(auto it = data.begin(); it != data.end(); it ++) {
        json t = *it;
        vector<double> p;
        p.push_back(t["px"]);
        p.push_back(t["py"]);
        points.push_back(Vector(p));
    }
}

json persistence_homology_split(json data)
{
	Points points;
	read_points_from_json(data, points);

	int max_d = 2;

    //Filtration* full_filtration = new RipsFiltration(points, max_d);
    //PersistentHomology sparse_rips_homology(full_filtration);
    Filtration* sparse_filtration = new SparseRipsFiltration(points, max_d, 1.0/3);
    PersistentHomology ph(sparse_filtration);

    std::vector<PHCycle> reduction;
    Filtration* filtration;
    int loaded_max_d;
    // TODO load these from saved data

    ph.compute_matrix(reduction);

	PersistenceDiagram *sparse_rips_pd = ph.compute_persistence(reduction, filtration, loaded_max_d);

	sparse_rips_pd->sort_pairs_by_persistence();

    std::stringstream result;
	for(unsigned i = 0; i < sparse_rips_pd->num_pairs(); i++)  {
		PersistentPair pairing = sparse_rips_pd->get_pair(i);
		//printf("%u %.7f %.7f\n", pairing.dim(), pairing.birth_time(), pairing.death_time());
        result << pairing.dim() << " "
            << pairing.birth_time() << " "
            << pairing.death_time() << "\n";
	}
    return result.str();
}

static void handle_query_call(struct mg_connection *c, struct http_message *hm) {

  json q = json::parse(string(hm->body.p, hm->body.len));
  json result = persistence_homology_split(q);

  /* Send result */
  std::string msg_content = result.dump();
  const std::string sep = "\r\n";

  std::stringstream ss;
  ss << "HTTP/1.1 200 OK"             << sep
    << "Content-Type: application/json" << sep
    << "Access-Control-Allow-Origin: *" << sep
    << "Content-Length: %d"             << sep << sep
    << "%s";

  mg_printf(c, ss.str().c_str(), (int) msg_content.size(), msg_content.c_str());
}

static void ev_handler(struct mg_connection *c, int ev, void *ev_data) {
  struct http_message *hm = (struct http_message *) ev_data;

  switch (ev) {
    case MG_EV_HTTP_REQUEST:
      if (mg_vcmp(&hm->uri, "/query") == 0) {
        handle_query_call(c, hm); /* Handle RESTful call */
      }
      else {
        mg_serve_http(c, hm, s_http_server_opts); /* Serve static content */
      }
      break;
    default:
      break;
  }
}

int main(int argc, char *argv[]) {

  struct mg_mgr mgr;
  struct mg_connection *c;

  mg_mgr_init(&mgr, NULL);
  c = mg_bind(&mgr, s_http_port, ev_handler);
  mg_set_protocol_http_websocket(c);

  s_http_server_opts.document_root = "../../webui";
  s_http_server_opts.enable_directory_listing = "yes";


  printf("Starting server on port %s\n", s_http_port);

  for (;;) {
    mg_mgr_poll(&mgr, 1000);
  }


  mg_mgr_free(&mgr);

  return 0;
}
