//
// Created by Dennis Goldfarb on 7/13/15.
//

#ifndef MSACQUISITIONSIMULATOR_PEPTIDEFORDIGESTION_H
#define MSACQUISITIONSIMULATOR_PEPTIDEFORDIGESTION_H


#include "Peptide.h"

class PeptideForDigestion : public Peptide {

private:

public:
	PeptideForDigestion() : Peptide() {}

	PeptideForDigestion(Protein* protein, int start, int end, std::map<int,PTM*> index2mod_for_pep,
						int n_cleavage_index, int c_cleavage_index, double abundance) :
	Peptide(protein, start, end, index2mod_for_pep), n_cleavage_index(n_cleavage_index),
	c_cleavage_index(c_cleavage_index) {
		this->abundance = abundance;
	}

	int n_cleavage_index;
	int c_cleavage_index;

	const char& operator[](std::size_t i) const;
};


#endif //MSACQUISITIONSIMULATOR_PEPTIDEFORDIGESTION_H
