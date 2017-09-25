#include "generalized_metric_space.h"

GeneralizedMetricSpace::GeneralizedMetricSpace(int _numPoints, double** _distances) : MetricSpace(_numPoints) {
	distances = _distances;
  this->numPoints = _numPoints;
}

GeneralizedMetricSpace::~GeneralizedMetricSpace()  {
	for(int i = 0; i < this->numPoints; i++)
		delete [] distances[i];
	delete [] distances;
}

double GeneralizedMetricSpace::distance(int _i, int _j)  {
	return distances[_i][_j];
}
