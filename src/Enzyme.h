//
// Created by Dennis Goldfarb on 6/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_ENZYME_H
#define MSACQUISITIONSIMULATOR_ENZYME_H

#include <iostream>
#include "PeptideForDigestion.h"

class Enzyme {
private:

public:
	Enzyme() {}
	Enzyme(std::string name) : name(name) {};

	std::string name;

	virtual double calc_cleavage_probability(const PeptideForDigestion &peptide);
	virtual double calc_cleavage_probability(const PeptideForDigestion &peptide, int index);

	Enzyme& operator=(const Enzyme& other) {
		return *this;
	};
};



#endif //MSACQUISITIONSIMULATOR_ENZYME_H
