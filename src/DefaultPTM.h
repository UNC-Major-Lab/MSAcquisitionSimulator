//
// Created by Dennis Goldfarb on 7/10/15.
//

#ifndef MSACQUISITIONSIMULATOR_DEFAULTPTM_H
#define MSACQUISITIONSIMULATOR_DEFAULTPTM_H

#include "PTM.h"
#include <math.h>

class DefaultPTM : public PTM {

private :

public:
	DefaultPTM() {}

	DefaultPTM(const std::string name, const std::string abbreviation, const std::string residues, const std::string formula, double site_probability,
			   double bind_energy, double bind_area,
			   double relative_abundance, bool does_block_cleavage, bool is_post_digestion) :
			PTM(name, abbreviation, residues, formula, bind_energy, bind_area, is_post_digestion),
			relative_abundance(relative_abundance), does_block_cleavage(does_block_cleavage),
			site_probability(site_probability) {

		log_relative_abundance = log(relative_abundance);
		log_one_minus_relative_abundance = log(1-relative_abundance);
	}

	DefaultPTM& operator=(const DefaultPTM& other) {
		return *this;
	};

	~DefaultPTM() {}

	double relative_abundance;
	double log_relative_abundance;
	double log_one_minus_relative_abundance;
	double site_probability;
	bool does_block_cleavage;

	std::vector<int> get_localizations(const std::string & sequence);
};


#endif //MSACQUISITIONSIMULATOR_DEFAULTPTM_H
