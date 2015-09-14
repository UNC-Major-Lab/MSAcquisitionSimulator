//
// Created by Dennis Goldfarb on 7/12/15.
//

#ifndef MSACQUISITIONSIMULATOR_PEPTIDE_H
#define MSACQUISITIONSIMULATOR_PEPTIDE_H

#include <vector>
#include <map>
#include <algorithm>
#include "Protein.h"
#include "PTM.h"
#include "Residue.h"

class Peptide {

private:

public:
	Peptide() {}

	Peptide(Protein* protein, int start, int end, std::map<int,PTM*> index2mod_for_pep) :
			protein(protein), start(start), end(end), index2mod_for_pep(index2mod_for_pep) {}

	Protein* protein;
	int start;
	int end;
	std::map<int,PTM*> index2mod_for_pep;
	double abundance=0;

	const char& operator[](std::size_t i);

	std::string get_sequence() const;
	std::string get_modified_sequence() const;
	int num_basic_residues() const;
	bool operator<(const Peptide& other) const;
	std::vector<unsigned int> get_composition() const;
	std::vector<unsigned int> get_composition(int charge) const;
	double calc_monoisotopic_mass();
	double monoisotopic_mass=0;

	static bool less_mz(const Peptide& a, const Peptide& b);
	bool operator==(const Peptide &other) const;

	static struct _CompareDoubleField
	{
		bool operator() (double left, const Peptide & right) {
			return left < right.monoisotopic_mass;
		}

	} CompareDoubleField;
};


#endif //MSACQUISITIONSIMULATOR_PEPTIDE_H
