#include "fixed_filtration.h"


FixedFiltration::FixedFiltration() : Filtration(2)  {
    Points points;

    for(int i = 0; i < 6; i ++) {
        vector<double> p;
        p.push_back(i);
        p.push_back(0);
        points.push_back(Vector(p));
    }

	this->num_points = points.size();
	this->distances = new double*[num_points];

	for(int i = 0; i < num_points; i++)  {
		distances[i] = new double[num_points];
		distances[i][i] = 0;
		for(int j = 0; j < i; j++)  {
			distances[i][j] = (points[i]-points[j]).l2Norm();
			distances[j][i] = distances[i][j];
		}
	}
	metric_space = new GeneralizedMetricSpace(num_points, distances);
}

FixedFiltration::~FixedFiltration()  {
	for(int i = 0; i < num_points; i++)
		delete [] distances[i];
	delete [] distances;
	delete metric_space;
}

void FixedFiltration::build_metric() {
	metric_space = new GeneralizedMetricSpace(num_points, distances);
}

bool FixedFiltration::build_filtration()  {
	if(all_simplices.size() > 0)
		all_simplices.clear();
	std::set<int> all_vertices;
	for(int i = 0; i < num_points; i++)
		all_vertices.insert(i);

    // all nodes are 0-simplices
    for(int i = 0; i < this->num_points; i ++) {
        std::vector<int> s = {i+1};
        this->all_simplices.push_back(Simplex(s, this->metric_space));
    }

    // add 1-simplices
    for(int i = 0; i < this->num_points; i ++) {
        for(int j = i+1; j < this->num_points; j ++) {
            if(distances[i][j] < 1.2) {
                std::vector<int> s = {i+1, j+1};
                this->all_simplices.push_back(Simplex(s, this->metric_space));
            }
        }
    }

    // sort by distance
	for(int s = 0; s < all_simplices.size(); s++)
		this->all_simplices[s].compute_simplex_distance();
	std::sort(this->all_simplices.begin(), this->all_simplices.end());

    // print simplex info
	std::vector<int> simplex_count(2,0);
	for(int i = 0; i < all_simplices.size(); i++)  {
		simplex_count[all_simplices[i].dim()]=simplex_count[all_simplices[i].dim()]+1;
	}
	for(int i = 0; i < simplex_count.size(); i++)  {
		long max_num_simplices = num_points;
		for(int j = 1; j <= i; j++)
			max_num_simplices *= (num_points-j);
		for(int j = 1; j <= i; j++)
			max_num_simplices /= (j+1);
		std::cout << "number of " << i << " simplices: " << simplex_count[i] << std::endl;
	}

	return true;
}
