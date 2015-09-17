//
// Created by Dennis Goldfarb on 6/29/15.
//


#include "MolecularFormula.h"


MolecularFormula::MolecularFormula(std::string formula) {
	boost::char_separator<char> sep(",");
	boost::tokenizer<boost::char_separator<char> > tokens(formula, sep);
	for (auto itr = tokens.begin(); itr != tokens.end(); itr++) {
		std::string elem = "";
		int count = 0;
		for (unsigned long i = 0; i < itr->size(); i++) {
			if (isalpha((*itr)[i])) {
				elem.push_back((*itr)[i]);
			} else {
				count = atoi((*itr).substr(i).c_str());
				break;
			}
		}
		if (elements::name2element.find(elem) != elements::name2element.end()) {
			element2count[elements::name2element.at(elem)] = count;
		} else if (subatomic_particles::name2particle.find(elem) != subatomic_particles::name2particle.end()) {
			particle2count[subatomic_particles::name2particle.at(elem)] = count;
		}
	}
}

double MolecularFormula::get_monoisotopic_mass() {
	double monoisotopic_mass = 0;
	for (std::pair<int, int> pair : element2count) {
		monoisotopic_mass += mercury::elemMasses[pair.first][0] * pair.second;
	}

	for (std::pair<const SubatomicParticle *, int> pair : particle2count) {
		monoisotopic_mass += pair.first->mass * pair.second;
	}

	return monoisotopic_mass;
}

double MolecularFormula::get_average_mass() {
	return 0;
}
