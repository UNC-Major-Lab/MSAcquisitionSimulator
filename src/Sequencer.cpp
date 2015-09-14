//
// Created by Dennis Goldfarb on 8/20/15.
//

#include "Sequencer.h"

void Sequencer::initialize(std::string fasta_in_path) {
	std::cout << "Parsing FASTA file and digesting proteins..." << std::endl;
	proteins = parse_FASTA(fasta_in_path.c_str());
	std::unordered_set<Peptide, KeyHasher> digested_peptides;

	for (Protein &protein : proteins) {
		digest_protein(&protein, digested_peptides);
	}

	for (auto peptide = digested_peptides.begin(); peptide != digested_peptides.end(); ++peptide) {
		peptides.push_back(*peptide);
	}

	std::sort(peptides.begin(), peptides.end(), Peptide::less_mz);
	std::cout << "Digestion complete." << std::endl;
}



void Sequencer::sequence_ms2_scan(MS2Scan *scan) {
	if (scan->peptide2intensity.size() == 0) return;
	double rand = g_uniform_distribution(g_rng)*std::max(scan->TIC, scan->target_total_ion_count);

	std::string selected_seq;
	double probability=0;
	double sum = 0;
	for (auto itr = scan->peptide2intensity.begin(); itr != scan->peptide2intensity.end(); ++itr) {
		sum+= itr->second;
		if (rand <= sum) {
			selected_seq = itr->first;
			break;
		}
	}

	std::set<std::string> candidate_peptides;
	for (int z = 2; z <= 4; z++) {
		double mass = (scan->precursor_peak.mz * z) - (z * mercury::elemMasses[(elements::ELEMENTS::H)][0]);
		std::vector<std::string> matched_peptides = get_peptides_within_range(mass - mass_tolerance,
																			  mass + mass_tolerance);
		for (std::string &pep : matched_peptides) candidate_peptides.insert(pep);
	}

	bool correct = false;
	if (candidate_peptides.find(selected_seq) != candidate_peptides.end()) {
		// correct sequence
		probability = scan->peptide2intensity[selected_seq] / std::max(scan->TIC, scan->target_total_ion_count);
		correct = true;
	} else {
		// incorrect sequence
		if (candidate_peptides.size() == 0 || g_uniform_distribution(g_rng) < .5) {
			selected_seq = "DECOY";
		} else {
			int random_index = (int) (std::floor(g_uniform_distribution(g_rng)*candidate_peptides.size()));
			random_index = std::min(random_index, (int) candidate_peptides.size()-1);

			auto itr = candidate_peptides.begin();
			for (int i = 0; i < random_index && itr != candidate_peptides.end(); i++) ++itr;

			selected_seq = *itr;
			probability = get_null_probability();
		}
	}
	scan->peptide = selected_seq;
	scan->probability = probability;

	if (pep2prot.find(selected_seq) != pep2prot.end()) {
		scan->proteins = pep2prot[selected_seq];
	}

	//if (correct)
	//	std::cout << scan->precursor_peak.mz << "\t" << selected_seq << "\t" << probability <<  "\t" << candidate_peptides.size() << "\t" << correct << std::endl;

}

void Sequencer::digest_protein(Protein *protein, std::unordered_set<Peptide, KeyHasher> &digested_peptides) {
	std::map<int, std::vector<DBPTM *>> index2mod;

	for (DBPTM* p : ptms) {
		std::vector<int> localizations = p->get_localizations(protein->sequence);
		for (int i : localizations) {
			index2mod[i].push_back(p);
		}
	}

	for (int start = 0; start < protein->sequence.size(); start++) {
		double mass = 0;
		int num_missed_cleavages = 0;
		int prev_num_enzymatic_termini = 2;

		for (DBEnzyme &e : enzymes) {
			switch (e.terminus) {
				case TERMINUS::N_TERM:
					if (e.residues.find(protein->sequence[start]) == e.residues.npos) {
						--prev_num_enzymatic_termini;
					}
					break;
				case TERMINUS::C_TERM:
					if (start > 0 && e.residues.find(protein->sequence[start-1]) == e.residues.npos) {
						--prev_num_enzymatic_termini;
					}
					break;
			}
		}
		if (prev_num_enzymatic_termini < min_enzymatic_termini) continue;

		for (int end = start; end < protein->sequence.size(); end++) {
			if (residues::name2residue.find(protein->sequence[end]) == residues::name2residue.end()) break;
			mass += residues::name2residue.at(protein->sequence[end])->monoisotopic_mass;

			if (mass < min_mass) continue;
			if (mass > max_mass) break;
			if ((end-start)+1 < min_aa) continue;

			int num_enzymatic_termini = prev_num_enzymatic_termini;

			for (DBEnzyme &e : enzymes) {
				switch (e.terminus) {
					case TERMINUS::N_TERM:
						if (e.residues.find(protein->sequence[end]) != e.residues.npos) {
							++num_missed_cleavages;
						}
						if (end < protein->sequence.size()-1 && e.residues.find(protein->sequence[end+1]) == e.residues.npos) {
							--num_enzymatic_termini;
						}
						break;
					case TERMINUS::C_TERM:
						if (e.residues.find(protein->sequence[end-1]) != e.residues.npos) {
							++num_missed_cleavages;
						}
						if (end < protein->sequence.size()-1 && e.residues.find(protein->sequence[end]) == e.residues.npos) {
							--num_enzymatic_termini;
						}
						break;
				}
			}


			if (num_missed_cleavages > max_missed_cleavages) break;
			if (num_enzymatic_termini < min_enzymatic_termini) continue;


			std::vector<std::vector<PTMLocation>> locations;
			for (int i = start; i < end; i++) {

				if (index2mod.find(i) != index2mod.end()) {
					locations.push_back(std::vector<PTMLocation>());
					bool is_static = false;
					for (DBPTM* ptm : index2mod[i]) {
						if (ptm->is_static) {
							is_static = true;
						}
						locations[locations.size()-1].push_back(PTMLocation(i, ptm));
					}
					if (!is_static)
						locations[locations.size()-1].push_back(PTMLocation(i, nullptr));
				}

			}

			Enumerator<PTMLocation> enumerator(locations, max_dynamic_mods);
			std::vector<short> state;

			while (enumerator.has_next()) {
				state = enumerator.next_combination();
				std::map<int, PTM*> peptide_index2mod;

				for (int i = 0; i < state.size(); i++) {
					PTMLocation pl = enumerator.states[i][state[i]];
					if (pl.ptm != nullptr) {
						peptide_index2mod[pl.index] = pl.ptm;
					}
				}

				Peptide p(protein, start, end, peptide_index2mod);
				p.monoisotopic_mass = p.calc_monoisotopic_mass();
				digested_peptides.insert(p);

				std::string mod_seq = p.get_modified_sequence();
				if (pep2prot.find(mod_seq) == pep2prot.end()) {
					pep2prot.insert(std::pair<std::string, std::vector<Protein *>>(mod_seq, std::vector<Protein *>()));
				}
				pep2prot[mod_seq].push_back(p.protein);
			}



		}
	}
}


double Sequencer::get_null_probability() {
	return std::min(std::max(null_probability_dist(g_rng),0.0),1.0);
}

std::vector<std::string> Sequencer::get_peptides_within_range(double low_mass, double high_mass) {
	std::vector<std::string> matched_peptides;

	auto low = std::upper_bound(peptides.begin(), peptides.end(), low_mass, Peptide::CompareDoubleField);

	for (; low != peptides.end() && low->monoisotopic_mass <= high_mass; ++low) {
		matched_peptides.push_back(low->get_modified_sequence());
	}

	return matched_peptides;
}
