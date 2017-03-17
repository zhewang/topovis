#ifndef COVER_H
#define COVER_H

#include "simplex.h"
#include "filtration.h"

#include<vector>
#include<map>


class Cover {
	public:
		Cover(Filtration *filtration, std::map<int,int> &labels);

        std::vector< std::vector<Simplex> > subComplexes;
        std::map<std::string, std::vector<Simplex> > blowupComplexes;
};

#endif
