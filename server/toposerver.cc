#include <iostream>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <typeinfo>
#include <cstdlib>
#include <fstream>

//#include <boost/random.hpp>
//#include <boost/generator_iterator.hpp>
//#include <boost/date_time/gregorian/gregorian.hpp>
//#include <boost/date_time/posix_time/posix_time.hpp>
//#include <boost/algorithm/string.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>

#include "mongoose.h"
#include "json.hpp"

#include "topology/persistence.h"
#include "geometry/covertree.h"
#include "topology/sparse_rips_filtration.h"
#include "topology/rips_filtration.h"
#include "topology/fixed_filtration.h"

using json = nlohmann::json;


static const char *s_http_port = "8800";
static struct mg_serve_http_opts s_http_server_opts;

void read_points_from_json(json& data, Points& points)
{
    for(auto it = data["points"].begin(); it != data["points"].end(); it ++) {
        json t = *it;
        vector<double> p;
        p.push_back(t["px"]);
        p.push_back(t["py"]);
        points.push_back(Vector(p));
    }
}

//json simple_ph()
//{
    //// build filtration
    //Filtration* f = new FixedFiltration();
    //f->build_filtration();

    //std::set<int> selection = {1,2,3};

    //Filtration* selectedComplex = new FixedFiltration();
    //Filtration* unselectedComplex = new FixedFiltration();
    //Filtration* blowupComplex = new FixedFiltration();

    //f->binarySplit(selection, selectedComplex, unselectedComplex, blowupComplex);
    //selectedComplex->print();
    //unselectedComplex->print();
    //blowupComplex->print();

    //// calculate ph
    //PersistentHomology ph(f);

    //std::vector<PHCycle> reduction;

    //ph.compute_matrix(reduction);

	//PersistenceDiagram *pd = ph.compute_persistence(reduction);

    //pd->sort_pairs_by_persistence();

    //std::stringstream result;
	//for(unsigned i = 0; i < pd->num_pairs(); i++)  {
		//PersistentPair pairing = pd->get_pair(i);
		////printf("%u %.7f %.7f\n", pairing.dim(), pairing.birth_time(), pairing.death_time());
        //result << pairing.dim() << " "
            //<< pairing.birth_time() << " "
            //<< pairing.death_time() << "\n";
	//}

    //delete pd;
    //return result.str();
//}

//json persistence_homology(json data)
//{
	//Points points;
	//read_points_from_json(data, points);

	//int max_d = 2;

    //Filtration* full_filtration = new RipsFiltration(points, max_d);
    //PersistentHomology ph(full_filtration);
    ////Filtration* sparse_filtration = new SparseRipsFiltration(points, max_d, 1.0/3);
    ////PersistentHomology ph(sparse_filtration);

    //std::vector<PHCycle> reduction;

    //string sel_id = data["sel_id"];
    //string filename = sel_id + "_reduction.txt";

    //// load serialized vector
    //{
        //std::ifstream ifs(filename);
        //boost::archive::text_iarchive ia(ifs);
        //ia & reduction;
    //}

	//PersistenceDiagram *rips_pd = ph.compute_persistence(reduction);

	//rips_pd->sort_pairs_by_persistence();

    //std::stringstream result;
	//for(unsigned i = 0; i < rips_pd->num_pairs(); i++)  {
		//PersistentPair pairing = rips_pd->get_pair(i);
		////printf("%u %.7f %.7f\n", pairing.dim(), pairing.birth_time(), pairing.death_time());
        //result << pairing.dim() << " "
            //<< pairing.birth_time() << " "
            //<< pairing.death_time() << "\n";
	//}

    //delete rips_pd;
    //return result.str();
//}

json compute_reduction_matrix(json data)
{
	Points points;
	read_points_from_json(data, points);

	int max_d = 2;

    //Filtration* full_filtration = new RipsFiltration(points, max_d);
    //PersistentHomology ph(full_filtration);
    Filtration* sparse_filtration = new SparseRipsFiltration(points, max_d, 1.0/3);
    sparse_filtration->build_filtration();
    std::vector<Simplex> sc = sparse_filtration->get_complex();

    std::vector<PHCycle> reduction = PersistentHomology::compute_matrix(sc);

    //string sel_id = data["sel_id"];
    //string filename = sel_id + "_reduction_basis.txt";

    //// serialize vector
    //{
        //std::ofstream ofs(filename);
        //boost::archive::text_oarchive oa(ofs);
        //oa & reduction;
    //}

    //std::vector<PHCycle> new_reduction;

    //// load serialized vector
    //{
        //std::ifstream ifs(filename);
        //boost::archive::text_iarchive ia(ifs);
        //ia & new_reduction;
    //}

	PersistenceDiagram *sparse_rips_pd = PersistentHomology::compute_persistence(reduction, sc);

	sparse_rips_pd->sort_pairs_by_persistence();

    std::stringstream result;
	for(unsigned i = 0; i < sparse_rips_pd->num_pairs(); i++)  {
		PersistentPair pairing = sparse_rips_pd->get_pair(i);
		//printf("%u %.7f %.7f\n", pairing.dim(), pairing.birth_time(), pairing.death_time());
        result << pairing.dim() << " "
            << pairing.birth_time() << " "
            << pairing.death_time() << "\n";
	}

    delete sparse_rips_pd;
    return result.str();
}

static void handle_query_call(struct mg_connection *c, struct http_message *hm) {

  json q = json::parse(string(hm->body.p, hm->body.len));
  json result = compute_reduction_matrix(q);
  //json result = persistence_homology(q);
  //json result = simple_ph();

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
