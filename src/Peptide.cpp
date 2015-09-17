//
// Created by Dennis Goldfarb on 7/12/15.
//

#include "Peptide.h"

const char& Peptide::operator[](std::size_t i) {
	return protein->sequence[i];
}

std::string Peptide::get_sequence() const {
	return protein->sequence.substr(start,(end-start)+1);
}


/* has fewer mods
 * if equal number of mods, has smaller lexicographical sequence
 * if equal sequences, has smaller lexicographical modified sequence
 */
bool Peptide::operator<(const Peptide &other) const {
	if (index2mod_for_pep.size() != other.index2mod_for_pep.size()) return index2mod_for_pep.size() < other.index2mod_for_pep.size();
	if (get_sequence() != other.get_sequence()) return get_sequence() < other.get_sequence();
 	return get_modified_sequence() < other.get_modified_sequence();
}

std::string Peptide::get_modified_sequence() const {
	std::string modified_sequence;
	int max = std::min(end, (int)protein->sequence.length()-1);
 	for (int i = start; i <= max; ++i) {
		if (index2mod_for_pep.find(i) != index2mod_for_pep.end()) {
			modified_sequence += index2mod_for_pep.at(i)->abbreviation;
		}
		modified_sequence.push_back(protein->sequence[i]);
	}
	return modified_sequence;
}

int Peptide::num_basic_residues() const {
	int basic_residues = 0;

	for (char c : get_sequence()) {
		if (c=='H' || c=='K' || c=='R') basic_residues++;
	}

	return basic_residues;
}

std::vector<unsigned int> Peptide::get_composition() const {
	std::vector<unsigned int> composition(mercury::MAX_ELEMENTS);

	for (char c : get_sequence()) {
		const Residue* r = residues::name2residue.at(c);
		for (std::pair<int, int> pair : r->molecular_formula.element2count) {
			composition[pair.first] += pair.second;
		}
	}

	for (std::pair<int,PTM*> pair : index2mod_for_pep) {
		for (std::pair<int, int> pair2 : pair.second->molecular_formula.element2count) {
			composition[pair2.first] += pair2.second;
		}
	}

	composition[elements::ELEMENTS::H] += 2;
	composition[elements::ELEMENTS::O] += 1;

	return composition;
}

std::vector<unsigned int> Peptide::get_composition(int charge) const {
	std::vector<unsigned int> composition = get_composition();
	composition[elements::ELEMENTS::H] += charge;
	return composition;
}

bool Peptide::less_mz(const Peptide &a, const Peptide &b) {
	return a.monoisotopic_mass < b.monoisotopic_mass;
}

double Peptide::calc_monoisotopic_mass() {
	if (monoisotopic_mass == 0) {
		std::vector<unsigned int> element_counts = get_composition();
		for (int i = 0; i < element_counts.size(); i++) {
			monoisotopic_mass += element_counts[i] * mercury::elemMasses[i][0];
		}
	}
	return monoisotopic_mass;
}

bool Peptide::operator==(const Peptide &other) const {
	return (get_modified_sequence() == other.get_modified_sequence());
}

