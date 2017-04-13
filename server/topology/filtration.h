#ifndef FILTRATION_H
#define FILTRATION_H

#include "../geometry/point_incs.h"
#include "simplex.h"
#include "simplicial_complex.h"

#include <set>
#include <vector>


class Filtration  {
	public:
		Filtration(int _maxD);
		~Filtration();

		Simplex get_simplex(int _t)  { return all_simplices[_t]; }
		int filtration_size()  { return all_simplices.size(); }
        SimplicialComplex get_complex();

		int maxD()  { return max_d; }

		virtual bool build_filtration();

        void addSimplex(Simplex s) { this->all_simplices.push_back(s); }

        void print();

	protected:
		std::vector<Simplex> all_simplices;

	private:
		int max_d;
};

#endif
