#include "simplex.h"

Simplex::Simplex()  {
	metric_space = NULL;
	cached_distance = 0;
    uid = "null";
}

Simplex::Simplex(const std::vector<int> & _simplex, MetricSpace* _metricSpace)  {
	simplex = _simplex;
	std::sort(simplex.begin(), simplex.end());

	metric_space = _metricSpace;
	cached_distance = 0;

    char unique_id[10*(simplex.size()+1)];
    sprintf(unique_id, "%u", simplex[0]);
    for(unsigned i = 1; i < simplex.size(); i++)
        sprintf(unique_id, "%s-%u", unique_id, simplex[i]);
    uid = std::string(unique_id);
}

Simplex::~Simplex()  {
}

std::vector<Simplex> Simplex::faces()  {
	std::vector<Simplex> all_faces;
	for(unsigned i = 0; i < simplex.size(); i++)  {
		std::vector<int> new_face;
		for(unsigned j = 0; j < simplex.size()-1; j++)  {
			int next_vertex = simplex[(j+i)%simplex.size()];
			new_face.push_back(next_vertex);
		}
		all_faces.push_back(Simplex(new_face,metric_space));
	}
	return all_faces;
}
