//
// Created by Dennis Goldfarb on 8/20/15.
//

#ifndef MSACQUISITIONSIMULATOR_SEQUENCER_H
#define MSACQUISITIONSIMULATOR_SEQUENCER_H


#include <string>
#include <unordered_set>
#include <map>
#include <iostream>
#include <vector>

#include "PTM.h"
#include "Enzyme.h"
#include "DBEnzyme.h"
#include "FASTAParser.h"
#include "PTMLocation.h"
#include "Enumerator.h"
#include "DBPTM.h"
#include "MS2Scan.h"
#include "Globals.h"

struct KeyHasher {
	std::size_t operator()(const Peptide& k) const
	{
		size_t my_hash = 0;
		my_hash = std::hash<std::string>()(k.get_modified_sequence());
		return my_hash;
	}
};

class Sequencer {
private:
	void initialize(std::string fasta_in_path);
	void digest_protein(Protein *protein, std::unordered_set<Peptide, KeyHasher> &digested_peptides);

	std::vector<std::string> get_peptides_within_range(double low_mass, double high_mass);
	double get_null_probability();
	boost::exponential_distribution<double> null_probability_dist;
public:
	Sequencer(std::string fasta_in_path, std::vector<DBPTM> ptms, std::vector<DBEnzyme> enzymes, double mass_tolerance,
			  int max_missed_cleavages, int min_enzymatic_termini, double min_mass, double max_mass, int max_dynamic_mods,
			  double null_lambda) :
			enzymes(enzymes), mass_tolerance(mass_tolerance),
			min_enzymatic_termini(min_enzymatic_termini), max_missed_cleavages(max_missed_cleavages),
			min_mass(min_mass), max_mass(max_mass), max_dynamic_mods(max_dynamic_mods), null_probability_dist(null_lambda) {
		for (DBPTM p : ptms) {
			this->ptms.push_back((new DBPTM(p)));
		}
		initialize(fasta_in_path);
	};

	~Sequencer() {
		for (DBPTM* p : ptms) delete p;
	}

	void sequence_ms2_scan(MS2Scan *scan);

	std::vector<DBPTM*> ptms;
	std::vector<DBEnzyme> enzymes;
	std::vector<Peptide> peptides;
	std::vector<Protein> proteins;
	std::map<std::string,std::vector<Protein*>> pep2prot;

	double mass_tolerance;
	int max_missed_cleavages;
	int min_enzymatic_termini;
	double min_mass;
	double max_mass;
	int max_dynamic_mods;
	int min_aa = 5;
};


#endif //MSACQUISITIONSIMULATOR_SEQUENCER_H
