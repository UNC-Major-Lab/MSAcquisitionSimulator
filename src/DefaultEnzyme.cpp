//
// Created by Dennis Goldfarb on 7/10/15.
//

#include "DefaultEnzyme.h"
#include "DefaultPTM.h"

double DefaultEnzyme::calc_log_no_cleavage_probability(const PeptideForDigestion &peptide, int index) {
	if (index == 0 || index == peptide.end) return 0;

	switch (terminus) {
		case TERMINUS::N_TERM:
			if (residues.find(peptide[index]) != residues.npos || residues[0] == '*') {
				if (blocking_residues.find(peptide[index-1]) == blocking_residues.npos) {
					if (peptide.index2mod_for_pep.find(index-1) == peptide.index2mod_for_pep.end()
						|| !((DefaultPTM*)peptide.index2mod_for_pep.at(index-1))->does_block_cleavage) {

						return log_no_cleavage_probability;
					}
				}
			}
			break;
		case TERMINUS::C_TERM:
			if (residues.find(peptide[index-1]) != residues.npos || residues[0] == '*') {
				if (blocking_residues.find(peptide[index]) == blocking_residues.npos) {
					if (peptide.index2mod_for_pep.find(index) == peptide.index2mod_for_pep.end()
						|| !((DefaultPTM*)peptide.index2mod_for_pep.at(index))->does_block_cleavage) {

						return log_no_cleavage_probability;
					}
				}
			}
			break;
	}

	return 0;
}

double DefaultEnzyme::calc_log_missed_cleavage_probability(const PeptideForDigestion &peptide, int num_missed) {
	switch (terminus) {
		case TERMINUS::N_TERM:
			for (const std::pair<int, PTM*> &p : peptide.index2mod_for_pep) {
				if (((DefaultPTM *) peptide.index2mod_for_pep.at(p.first))->does_block_cleavage) {
					if (residues.find(peptide[p.first]) != residues.npos || residues[0] == '*') {
						num_missed--;
					}
				}
			}
			break;
		case TERMINUS::C_TERM:
			for (const std::pair<int, PTM*> &p : peptide.index2mod_for_pep) {
				if (((DefaultPTM *) peptide.index2mod_for_pep.at(p.first))->does_block_cleavage) {
					if (residues.find(peptide[p.first - 1]) != residues.npos || residues[0] == '*') {
						num_missed--;
					}
				}
			}
			break;
	}
	return log_no_cleavage_probability*num_missed;
}

