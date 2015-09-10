//
// Created by Dennis Goldfarb on 7/10/15.
//

#ifndef MSACQUISITIONSIMULATOR_DEFAULTENZYME_H
#define MSACQUISITIONSIMULATOR_DEFAULTENZYME_H

#include "Enzyme.h"
#include "Globals.h"
#include "DefaultPTM.h"

class DefaultEnzyme : public Enzyme {

private:

public:
	DefaultEnzyme() {}
	DefaultEnzyme(std::string name, std::string residues, std::string blocking_residues, TERMINUS terminus, double cleavage_probability) :
	Enzyme(name), residues(residues), blocking_residues(blocking_residues), terminus(terminus), cleavage_probability(cleavage_probability) {
		log_cleavage_probability = log(cleavage_probability);
		log_no_cleavage_probability = log(1-cleavage_probability);
	};

	std::string residues;
	std::string blocking_residues;
	TERMINUS terminus;
	double cleavage_probability;
	double log_cleavage_probability;
	double log_no_cleavage_probability;

	double calc_log_missed_cleavage_probability(const PeptideForDigestion &peptide, int num_missed);
	double calc_log_no_cleavage_probability(const PeptideForDigestion &peptide, int index);


	DefaultEnzyme& operator=(const DefaultEnzyme& other) {
		return *this;
	};
};


#endif //MSACQUISITIONSIMULATOR_DEFAULTENZYME_H
