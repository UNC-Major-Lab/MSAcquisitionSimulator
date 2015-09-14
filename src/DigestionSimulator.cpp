//
// Created by Dennis Goldfarb on 7/12/15.
//

#include <unordered_set>
#include "DigestionSimulator.h"
#include "DefaultEnzyme.h"

double DigestionSimulator::digest_peptides(vector<PeptideForDigestion> &peptides, vector<int> &num_missed, vector<Peptide> &digested_peptides) {

	double min_prob_no_cleavage = 0;
	map<map<int,PTM*>, PeptideDigested> tmp_peptides;

	for (PeptideForDigestion &pep : peptides) {
		PeptideDigested dp = digest_peptide(pep, num_missed);
		min_prob_no_cleavage = min(min_prob_no_cleavage, dp.log_no_cleavage_prob);
		if (dp.abundance < PRUNE_THRESHOLD) continue;
		if (tmp_peptides.find(dp.index2mod_for_pep) == tmp_peptides.end()) tmp_peptides[dp.index2mod_for_pep] = dp;
		else tmp_peptides.at(dp.index2mod_for_pep).abundance += dp.abundance;

	}

	for (const pair<map<int,PTM*>, PeptideDigested> &p : tmp_peptides) {
		digested_peptides.push_back(p.second);
	}

	min_prob_no_cleavage = exp(min_prob_no_cleavage);

	return min_prob_no_cleavage;
}

PeptideDigested DigestionSimulator::digest_peptide(PeptideForDigestion &peptide, vector<int> &num_missed) {
	PeptideDigested digested_peptide(peptide);
	for (int i = 0; i < enzymes.size(); i++) {
		DefaultEnzyme* e = (DefaultEnzyme *) enzymes[i];
		digested_peptide.log_N_cleavage_prob += e->calc_log_no_cleavage_probability(peptide,peptide.n_cleavage_index);
		digested_peptide.log_C_cleavage_prob += e->calc_log_no_cleavage_probability(peptide,peptide.c_cleavage_index);
		digested_peptide.log_no_cleavage_prob += e->calc_log_missed_cleavage_probability(peptide, num_missed[i]);
	}

	if (peptide.n_cleavage_index != 0) {
		digested_peptide.log_N_cleavage_prob = log(1 - exp(digested_peptide.log_N_cleavage_prob));
	}
	if (peptide.c_cleavage_index != peptide.end) {
		digested_peptide.log_C_cleavage_prob = log(1 - exp(digested_peptide.log_C_cleavage_prob));
	}

	digested_peptide.abundance *= digested_peptide.calc_digestion_probability();

	return digested_peptide;
}
