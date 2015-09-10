//
// Created by Dennis Goldfarb on 7/12/15.
//

#ifndef MSACQUISITIONSIMULATOR_MODIFICATIONSIMULATOR_H
#define MSACQUISITIONSIMULATOR_MODIFICATIONSIMULATOR_H


#include "PTM.h"
#include "Peptide.h"
#include "Protein.h"
#include "PeptideForDigestion.h"
#include "PTMLocation.h"
#include <vector>
#include <map>
#include "Globals.h"

class ModificationSimulator {

private:

	void append_probabilities(std::vector<std::vector<std::pair<double,PTMLocation>>> &probabilities, std::map<int, std::vector<PTM *>> &index2mod, int i);
	static bool comp(const std::pair<double,PTMLocation>& a, const std::pair<double,PTMLocation>& b);

public:

	std::vector<PTM*> pre_digestion_modifications;
	std::vector<PTM*> post_digestion_modifications;

	ModificationSimulator() {}

	ModificationSimulator(std::vector<PTM*> pre_digestion_modifications, std::vector<PTM*> post_digestion_modifications) :
	pre_digestion_modifications(pre_digestion_modifications), post_digestion_modifications(post_digestion_modifications) {}

	~ModificationSimulator() {
		for (int i = 0; i < pre_digestion_modifications.size(); ++i) delete pre_digestion_modifications[i];
		for (int i = 0; i < post_digestion_modifications.size(); ++i) delete post_digestion_modifications[i];
	}

	std::map<int, std::vector<PTM*> > get_pre_digestion_localizations(const std::string& sequence);
	std::vector<Peptide> add_post_digestion_modifications(Peptide& p);

	std::vector<PeptideForDigestion> get_modified_peptides_for_digestion(Protein& p, std::map<int,std::vector<PTM*>>& index2mod,
																		 std::vector<std::vector<std::pair<double,PTMLocation>>> &probabilities,
																		 int start, int end);
};


#endif //MSACQUISITIONSIMULATOR_MODIFICATIONSIMULATOR_H
