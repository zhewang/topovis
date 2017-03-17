#ifndef FILTRATION_H
#define FILTRATION_H

#include "../geometry/point_incs.h"
#include "simplex.h"

#include <set>
#include <vector>

typedef std::vector<Simplex> Complex;

class Filtration  {
	public:
		Filtration(int _maxD);
		~Filtration();

		Simplex get_simplex(int _t)  { return all_simplices[_t]; }
		int filtration_size()  { return all_simplices.size(); }
        Complex get_complex() { return all_simplices; };

		int maxD()  { return max_d; }

		virtual bool build_filtration();

        bool binarySplit(std::set<int> &selectedVertices,
                       Filtration *selectedComplex,
                       Filtration *unselectedComplex,
                       Filtration *blowupComplex);

        void addSimplex(Simplex s) { this->all_simplices.push_back(s); }

        void print();

	protected:
		std::vector<Simplex> all_simplices;

	private:
		int max_d;
};

#endif
