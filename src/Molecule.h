//
// Created by Dennis Goldfarb on 7/10/15.
//

#ifndef MSACQUISITIONSIMULATOR_MOLECULE_H
#define MSACQUISITIONSIMULATOR_MOLECULE_H


#include "MolecularFormula.h"

class Molecule {
private:

public:
	MolecularFormula molecular_formula;

	Molecule(std::string formula) : molecular_formula(formula) {};

	Molecule() {}
};


#endif //MSACQUISITIONSIMULATOR_MOLECULE_H
