//
// Created by Dennis Goldfarb on 7/12/15.
//

#ifndef MSACQUISITIONSIMULATOR_PEPTIDEGENERATOR_H
#define MSACQUISITIONSIMULATOR_PEPTIDEGENERATOR_H

#include <vector>
#include "ModificationSimulator.h"
#include "DigestionSimulator.h"
#include "Peptide.h"
#include "Protein.h"
#include "DefaultEnzyme.h"
#include "DefaultPTM.h"
#include "Globals.h"

using namespace std;

class PeptideGenerator {

private:


public:
	PeptideGenerator() {}

	PeptideGenerator(std::vector<DefaultPTM>& ptms, std::vector<DefaultEnzyme>& enzymes);

	PeptideGenerator(ModificationSimulator modification_simulator, DigestionSimulator digestion_simulator) :
			modification_simulator(modification_simulator), digestion_simulator(digestion_simulator) {}

	vector<Peptide> generate_peptides(Protein &p);

	void add_enzyme(Enzyme* enzyme);
	void add_modification(PTM* ptm);

	ModificationSimulator modification_simulator;
	DigestionSimulator digestion_simulator;
};


#endif //MSACQUISITIONSIMULATOR_PEPTIDEGENERATOR_H
