//
// Created by Dennis Goldfarb on 7/22/15.
//

#ifndef MSACQUISITIONSIMULATOR_PTMLOCATION_H
#define MSACQUISITIONSIMULATOR_PTMLOCATION_H


#include "PTM.h"

struct PTMLocation {
	int index;
	PTM* ptm;

	PTMLocation(int index, PTM* ptm) : index(index), ptm(ptm) {}
	~PTMLocation() {}
};


#endif //MSACQUISITIONSIMULATOR_PTMLOCATION_H
