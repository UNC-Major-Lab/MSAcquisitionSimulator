//
// Created by Dennis Goldfarb on 7/12/15.
//

#include "PeptideGenerator.h"
#include "Residue.h"

vector<Peptide> PeptideGenerator::generate_peptides(Protein &p) {
	map<int, vector<PTM*>> index2mod = modification_simulator.get_pre_digestion_localizations(p.sequence);

	vector<Peptide> peptides;
	vector<vector<pair<double,PTMLocation>>> probabilities;
	vector<int> num_missed(digestion_simulator.enzymes.size(), 0);
	vector<PeptideForDigestion> peptides_for_digestion;
	double current_mass;
	int num_aa;

	for (int start = 0; start < p.sequence.length(); ++start) {

		for (int i = 0; i < num_missed.size(); ++i) num_missed[i] = 0;
		current_mass = 0;
		num_aa = 0;

		for (int end = start; end < p.sequence.length(); ++end) {
			++num_aa;
			if (residues::name2residue.find(p.sequence[end]) == residues::name2residue.end()) break;
			current_mass += residues::name2residue.at(p.sequence[end])->monoisotopic_mass;
			if (current_mass > MAX_MASS) break;

			peptides_for_digestion = modification_simulator.get_modified_peptides_for_digestion(p, index2mod,
																								probabilities,
																								start, end);
			if (num_aa > 1) {
				for (int i = 0; i < digestion_simulator.enzymes.size(); i++) {
					DefaultEnzyme* e = (DefaultEnzyme *) digestion_simulator.enzymes[i];
					switch (e->terminus) {
						case TERMINUS::N_TERM:
							if (e->residues.find(p.sequence[end]) != e->residues.npos || e->residues[0] == '*') {
								++num_missed[i];
							}
							break;
						case TERMINUS::C_TERM:
							if (e->residues.find(p.sequence[end-1]) != e->residues.npos || e->residues[0] == '*') {
								++num_missed[i];
							}
							break;
					}

				}
			}

			double min_prob_no_cleavage = digestion_simulator.digest_peptides(peptides_for_digestion, num_missed, peptides);
			if (min_prob_no_cleavage * p.abundance < PRUNE_THRESHOLD) break;

		}
	}

	/*double max_abundance = 0;
	double min_abundance = INFINITY;

	for (Peptide &p : peptides) {
		if (p.abundance > max_abundance) max_abundance = p.abundance;
		if (p.abundance < min_abundance) min_abundance = p.abundance;
	}

	cout << endl << max_abundance << " " << min_abundance << endl;*/

	return peptides;
}

void PeptideGenerator::add_modification(PTM* ptm) {
	if (ptm->is_post_digestion) modification_simulator.post_digestion_modifications.push_back(ptm);
	else modification_simulator.pre_digestion_modifications.push_back(ptm);
}

void PeptideGenerator::add_enzyme(Enzyme* enzyme) {
	digestion_simulator.enzymes.push_back(enzyme);
}

PeptideGenerator::PeptideGenerator(vector<DefaultPTM>& ptms, vector<DefaultEnzyme>& enzymes) {
	for (DefaultPTM ptm : ptms) add_modification(new DefaultPTM(ptm));
	for (DefaultEnzyme enzyme : enzymes) add_enzyme(new DefaultEnzyme(enzyme));
}
