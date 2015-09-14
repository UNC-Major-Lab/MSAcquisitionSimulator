//
// Created by Dennis Goldfarb on 7/23/15.
//

#ifndef MSACQUISITIONSIMULATOR_PEPTIDEDIGESTED_H
#define MSACQUISITIONSIMULATOR_PEPTIDEDIGESTED_H


#include "Peptide.h"
#include "PeptideForDigestion.h"

class PeptideDigested : public Peptide {

	public:

	PeptideDigested() : Peptide() {}

	PeptideDigested(PeptideForDigestion peptide) : Peptide(peptide) {
		start = peptide.n_cleavage_index;
		end = peptide.c_cleavage_index-1;
		abundance = peptide.abundance;
		protein = peptide.protein;
		index2mod_for_pep = peptide.index2mod_for_pep;

		for (int i = peptide.start; i < start; i++) {
			if (index2mod_for_pep.find(i) != index2mod_for_pep.end()) index2mod_for_pep.erase(i);
		}
		for (int i = end+1; i <= peptide.end; i++) {
			if (index2mod_for_pep.find(i) != index2mod_for_pep.end()) index2mod_for_pep.erase(i);
		}
	};

	double log_N_cleavage_prob = 0;
	double log_C_cleavage_prob = 0;
	double log_no_cleavage_prob = 0;

	double calc_digestion_probability();


};


#endif //MSACQUISITIONSIMULATOR_PEPTIDEDIGESTED_H
