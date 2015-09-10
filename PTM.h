//
// Created by Dennis Goldfarb on 6/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_PTM_H
#define MSACQUISITIONSIMULATOR_PTM_H


#include <iostream>
#include <vector>
#include <set>
#include "Molecule.h"


class PTM : public Molecule {

private:

public :

	PTM() {}

	PTM(const std::string name, const std::string abbreviation, const std::string residues, const std::string formula) :
			name(name), residues(residues), Molecule(formula), abbreviation(abbreviation), bind_energy(0), bind_area(0), is_post_digestion(false) {

	}

	PTM(const std::string name, const std::string abbreviation, const std::string residues, const std::string formula, double bind_energy, double bind_area, bool is_post_digestion) :
			name(name), residues(residues), bind_energy(bind_energy), bind_area(bind_area), is_post_digestion(is_post_digestion), Molecule(formula), abbreviation(abbreviation) {

	}

	PTM& operator=(const PTM& other) {
		return *this;
	};

	virtual ~PTM() {}

	double bind_energy;
	double bind_area;
	bool is_post_digestion;
	const std::string name;
	const std::string residues;
	const std::string abbreviation;

	virtual std::vector<int> get_localizations(const std::string & sequence);
};



#endif //MSACQUISITIONSIMULATOR_PTM_H
