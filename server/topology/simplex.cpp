#include "simplex.h"

Simplex::Simplex()  {
	metric_space = NULL;
	cached_distance = 0;
    uid = "null";
}

Simplex::Simplex(const std::vector<int> & _simplex, MetricSpace* _metricSpace)  {
	simplex = _simplex;
	std::sort(simplex.begin(), simplex.end());

  for(auto e : _simplex) this->loc_map[e] = e;

	metric_space = _metricSpace;
  cached_distance = -1;
  //compute_simplex_distance();
  compute_uid();
}

Simplex::Simplex(const std::vector<int> & _simplex, Points &points, std::map<int,int> _loc_map)  {
	simplex = _simplex;
	std::sort(simplex.begin(), simplex.end());

  this->loc_map = _loc_map;

	metric_space = NULL;
  cached_distance = -1;
  //compute_simplex_distance_from_points(points);
  compute_uid();
}

Simplex::Simplex(const std::vector<int> & _simplex, double _distance, std::map<int,int> _loc_map)  {
  simplex = _simplex;
  std::sort(simplex.begin(), simplex.end());

  this->loc_map = _loc_map;

  metric_space = NULL;
  cached_distance = _distance;

  compute_uid();
}

Simplex::~Simplex()  {
}

std::vector<Simplex> Simplex::faces(double** distances)  {
	std::vector<Simplex> all_faces;
	for(unsigned i = 0; i < simplex.size(); i++)  {
		std::vector<int> new_face;
		for(unsigned j = 0; j < simplex.size()-1; j++)  {
			int next_vertex = simplex[(j+i)%simplex.size()];
			new_face.push_back(next_vertex);
		}

    std::map<int,int> new_loc_map;
    for(auto e : new_face) new_loc_map[e] = this->loc_map[e];

    if(distances == NULL) {
      all_faces.push_back(Simplex(new_face, metric_space));
    } else {
      double new_distance = 0;
      for(unsigned i = 0; i < new_face.size(); i++)  {
        for(unsigned j = 0; j < i; j++)  {
          int loc_i = this->loc_map[new_face[i]];
          int loc_j = this->loc_map[new_face[j]];
          double next_dist = distances[loc_i][loc_j];
          new_distance = next_dist > new_distance ? next_dist : new_distance;
        }
      }
      all_faces.push_back(Simplex(new_face, new_distance, new_loc_map));
    }
	}
	return all_faces;
}

std::map<std::string, int> global_compare::order_map = std::map<std::string, int>();
