#ifndef COVER_H
#define COVER_H

#include "simplex.h"
#include "filtration.h"

#include<vector>
#include<map>


class Cover {
	public:
		Cover(Filtration *filtration, std::map<int,int> &labels);

        std::vector<Complex> subComplexes;
        std::map<std::string, Complex> blowupComplexes;
};

#endif
