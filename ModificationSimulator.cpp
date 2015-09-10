//
// Created by Dennis Goldfarb on 7/12/15.
//

#include "ModificationSimulator.h"
#include "DigestionSimulator.h"
#include "DefaultPTM.h"
#include "JointMultinomialEnumerator.h"
#include <math.h>

std::map<int, std::vector<PTM*> > ModificationSimulator::get_pre_digestion_localizations(const std::string &sequence) {
	std::map<int,std::vector<PTM*> > index2mod;
	for (PTM* p : pre_digestion_modifications) {
		std::vector<int> localizations = p->get_localizations(sequence);
		for (int i : localizations) {
			index2mod[i].push_back(p);
		}
	}
	return index2mod;
}

std::vector<PeptideForDigestion> ModificationSimulator::get_modified_peptides_for_digestion(Protein &p,
																							std::map<int, std::vector<PTM *>> &index2mod,
																							vector<vector<pair<double,PTMLocation>>> &probabilities,
																							int start, int end) {

	int start_extended = std::max(0, start-DigestionSimulator::DIGESTION_DISTANCE); // extend to the left at most 4 amino acids for determining probability of cleavage
	int end_extended = std::min(end+DigestionSimulator::DIGESTION_DISTANCE, (int) p.sequence.length()-1); // extend to the right at most 4 amino acids for determining probability of cleavage

	// create each permutation of modified peptides from start_extended to end_extended
	if (start == end) { 										// need to initialize probabilities
		probabilities.clear();
		for (int i = start_extended; i <= end_extended; ++i) {
			if (index2mod.find(i) != index2mod.end()) {
				append_probabilities(probabilities, index2mod, i);
			}
		}
	} else if (index2mod.find(end_extended) != index2mod.end() && end+DigestionSimulator::DIGESTION_DISTANCE < p.sequence.length()) {
		append_probabilities(probabilities, index2mod, end_extended);
	}
	std::vector<PeptideForDigestion> peptides;
	if (probabilities.size() > 0) {
		JointMultinomialEnumerator<PTMLocation> jme(probabilities, PRUNE_THRESHOLD / p.abundance);

		while (jme.has_next()) {
			pair<vector<short>, double> state = jme.next_combination();
			std::map<int,PTM*> index2mod_for_pep;
			for (int i = 0; i < state.first.size(); i++) {
				PTMLocation pl = jme.log_probabilities[i][state.first[i]].second;
				if (pl.ptm != nullptr) index2mod_for_pep[pl.index] = pl.ptm;
			}
			peptides.push_back(PeptideForDigestion(&p, start_extended, end_extended, index2mod_for_pep, start, end+1, p.abundance * state.second));
		}
	} else {
		peptides.push_back(PeptideForDigestion(&p, start_extended, end_extended, std::map<int,PTM*>(), start, end+1, p.abundance));
	}

	return peptides;
}

void ModificationSimulator::append_probabilities(std::vector<std::vector<std::pair<double, PTMLocation>>> &probabilities,
												 std::map<int, std::vector<PTM *>> &index2mod, int i) {

	probabilities.push_back(std::vector<pair<double,PTMLocation>>());
	vector<pair<double,PTMLocation>> &v = probabilities[probabilities.size()-1];
	double total_prob = 0;
	for (PTM* mod : index2mod[i]) {
		double relative_abundance = ((DefaultPTM*) mod)->relative_abundance;
		total_prob += relative_abundance;
		v.push_back({relative_abundance, {i,mod}});
	}
	if (total_prob < 1) {                        	// add no modification option
		v.push_back({1-total_prob, {i,nullptr}});
	} else { 										// normalize to sum = 1
		for (pair<double,PTMLocation> &prob : v) {
			prob.first /= total_prob;
		}
	}

	sort(v.begin(), v.end(), comp);
	for (pair<double, PTMLocation> &p : v) {
		p.first = log(p.first);
	}
}

bool ModificationSimulator::comp(const pair<double,PTMLocation> &a, const pair<double,PTMLocation> &b) {
	return a.first > b.first;
}