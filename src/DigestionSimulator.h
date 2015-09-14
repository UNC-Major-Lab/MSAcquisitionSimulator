//
// Created by Dennis Goldfarb on 7/12/15.
//

#ifndef MSACQUISITIONSIMULATOR_DIGESTIONSIMULATOR_H
#define MSACQUISITIONSIMULATOR_DIGESTIONSIMULATOR_H


#include <vector>
#include "Enzyme.h"
#include "Peptide.h"
#include "PeptideForDigestion.h"
#include "PeptideDigested.h"

using namespace std;

class DigestionSimulator {

private:
	PeptideDigested digest_peptide(PeptideForDigestion &peptide, vector<int> &num_missed);

public:

	static const int DIGESTION_DISTANCE = 4;

	vector<Enzyme*> enzymes;

	DigestionSimulator() {}
	DigestionSimulator(vector<Enzyme*> enzymes) : enzymes(enzymes) {}

	double digest_peptides(vector<PeptideForDigestion> & peptides, vector<int> &num_missed, vector<Peptide> &digested_peptides);
};


#endif //MSACQUISITIONSIMULATOR_DIGESTIONSIMULATOR_H
