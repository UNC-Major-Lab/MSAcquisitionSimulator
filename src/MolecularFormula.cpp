//
// Created by Dennis Goldfarb on 6/29/15.
//


#include "MolecularFormula.h"


MolecularFormula::MolecularFormula(std::string formula) {
	std::map<const std::string, int> name2element {{"H",  elements::ELEMENTS::H},
												   {"C",  elements::ELEMENTS::C},
												   {"N",  elements::ELEMENTS::N},
												   {"O",  elements::ELEMENTS::O},
												   {"S",  elements::ELEMENTS::S},
												   {"Br", elements::ELEMENTS::Br},
												   {"Cl", elements::ELEMENTS::Cl},
												   {"I",  elements::ELEMENTS::I},
												   {"P",  elements::ELEMENTS::P},
												   {"F",  elements::ELEMENTS::F},
												   {"Si", elements::ELEMENTS::Si},
												   {"Fe", elements::ELEMENTS::Fe},
												   {"B",  elements::ELEMENTS::B},
												   {"Se", elements::ELEMENTS::Se},
												   {"Co", elements::ELEMENTS::Co},
												   {"Mg", elements::ELEMENTS::Mg},
												   {"Au", elements::ELEMENTS::Au},
												   {"As", elements::ELEMENTS::As},
												   {"Sn", elements::ELEMENTS::Sn},
												   {"Na", elements::ELEMENTS::Na},
												   {"K",  elements::ELEMENTS::K},
												   {"Te", elements::ELEMENTS::Te},
												   {"Zn", elements::ELEMENTS::Zn},
												   {"Ge", elements::ELEMENTS::Ge},
												   {"Ca", elements::ELEMENTS::Ca},
												   {"Sb", elements::ELEMENTS::Sb},
												   {"Cu", elements::ELEMENTS::Cu},
												   {"Al", elements::ELEMENTS::Al},
												   {"Mn", elements::ELEMENTS::Mn},
												   {"Pt", elements::ELEMENTS::Pt},
												   {"Gd", elements::ELEMENTS::Gd},
												   {"Hg", elements::ELEMENTS::Hg},
												   {"Mo", elements::ELEMENTS::Mo},
												   {"Sr", elements::ELEMENTS::Sr},
												   {"Ga", elements::ELEMENTS::Ga},
												   {"Ni", elements::ELEMENTS::Ni},
												   {"Pb", elements::ELEMENTS::Pb},
												   {"Ag", elements::ELEMENTS::Ag},
												   {"Bi", elements::ELEMENTS::Bi},
												   {"Tl", elements::ELEMENTS::Tl},
												   {"Cr", elements::ELEMENTS::Cr},
												   {"Rb", elements::ELEMENTS::Rb},
												   {"Zr", elements::ELEMENTS::Zr},
												   {"Ti", elements::ELEMENTS::Ti},
												   {"W",  elements::ELEMENTS::W},
												   {"Be", elements::ELEMENTS::Be},
												   {"V",  elements::ELEMENTS::V},
												   {"Cd", elements::ELEMENTS::Cd},
												   {"Ba", elements::ELEMENTS::Ba},
												   {"Ta", elements::ELEMENTS::Ta},
												   {"Li", elements::ELEMENTS::Li},
												   {"Cs", elements::ELEMENTS::Cs},
												   {"Pd", elements::ELEMENTS::Pd},
												   {"Ce", elements::ELEMENTS::Ce},
												   {"Ru", elements::ELEMENTS::Ru},
												   {"La", elements::ELEMENTS::La},
												   {"Nd", elements::ELEMENTS::Nd},
												   {"Re", elements::ELEMENTS::Re},
												   {"Hf", elements::ELEMENTS::Hf},
												   {"Th", elements::ELEMENTS::Th},
												   {"He", elements::ELEMENTS::He},
												   {"Ar", elements::ELEMENTS::Ar},
												   {"Lu", elements::ELEMENTS::Lu},
												   {"U",  elements::ELEMENTS::U},
												   {"Kr", elements::ELEMENTS::Kr},
												   {"Ir", elements::ELEMENTS::Ir},
												   {"In", elements::ELEMENTS::In},
												   {"Rh", elements::ELEMENTS::Rh},
												   {"Ho", elements::ELEMENTS::Ho},
												   {"Dy", elements::ELEMENTS::Dy},
												   {"Yb", elements::ELEMENTS::Yb},
												   {"Eu", elements::ELEMENTS::Eu},
												   {"Os", elements::ELEMENTS::Os},
												   {"Pr", elements::ELEMENTS::Pr},
												   {"Tb", elements::ELEMENTS::Tb},
												   {"Er", elements::ELEMENTS::Er},
												   {"Xe", elements::ELEMENTS::Xe},
												   {"Sc", elements::ELEMENTS::Sc},
												   {"Ne", elements::ELEMENTS::Ne},
												   {"Sm", elements::ELEMENTS::Sm},
												   {"Tm", elements::ELEMENTS::Tm},
												   {"Nb", elements::ELEMENTS::Nb}};

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
		if (name2element.find(elem) != name2element.end()) {
			element2count[name2element.at(elem)] = count;
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
