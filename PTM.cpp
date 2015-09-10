//
// Created by Dennis Goldfarb on 6/29/15.
//

#include "PTM.h"

std::vector<int> PTM::get_localizations(const std::string &sequence) {
	std::vector<int> indices;

	// N-terminal mods
	if (residues.find('n') != residues.npos) indices.push_back(0);

	// Residue mods
	for (int i = 0; i < sequence.length(); i++) {
		if (residues.find(sequence[i]) != residues.npos) indices.push_back(i);
	}

	// C-terminal mods
	if (residues.find('c') != residues.npos) indices.push_back((int) sequence.length()-1);

	return indices;
}
