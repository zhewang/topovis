#ifndef FIXEDFILTRATION_H
#define FIXEDFILTRATION_H

#include "filtration.h"
#include "generalized_metric_space.h"

class FixedFiltration : public Filtration  {
	public:
		FixedFiltration();
		~FixedFiltration();

		virtual bool build_filtration();

        void build_metric();

	protected:
		int num_points;
		double** distances;
		MetricSpace* metric_space;

	private:
};

#endif
