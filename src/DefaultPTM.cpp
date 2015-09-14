//
// Created by Dennis Goldfarb on 7/10/15.
//

#include "DefaultPTM.h"
#include "Globals.h"

std::vector<int> DefaultPTM::get_localizations(const std::string & sequence)  {
	std::vector<int> indices;

	// N-terminal mods
	if (residues.find('n') != residues.npos) {
		if (g_uniform_distribution(g_rng) <= site_probability) {
			indices.push_back(0);
		}
	}

	// Residue mods
	for (int i = 0; i < sequence.length(); i++) {
		if (residues.find(sequence[i]) != residues.npos) {
			if (g_uniform_distribution(g_rng) <= site_probability) {
				indices.push_back(i);
			}
		}
	}

	// C-terminal mods
	if (residues.find('c') != residues.npos) {
		if (g_uniform_distribution(g_rng) <= site_probability) {
			indices.push_back((int) sequence.length()-1);
		}
	}
	return indices;
}