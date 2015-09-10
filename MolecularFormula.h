//
// Created by Dennis Goldfarb on 6/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_MOLECULARFORMULA_H
#define MSACQUISITIONSIMULATOR_MOLECULARFORMULA_H


#include <iostream>
#include <map>
#include "Element.h"
#include "SubatomicParticle.h"
#include "libmercury++.h"

class MolecularFormula {
public:

	MolecularFormula() {

	}

	MolecularFormula(std::string formula);

	MolecularFormula(std::map<const elements::ELEMENTS, int> element2count, std::map<const SubatomicParticle *, int> particle2count) :
			element2count(element2count), particle2count(particle2count) {

	}

	double get_monoisotopic_mass();
	double get_average_mass();

	std::map<const elements::ELEMENTS, int> element2count;
	std::map<const SubatomicParticle *, int> particle2count;
};



#endif //MSACQUISITIONSIMULATOR_MOLECULARFORMULA_H
