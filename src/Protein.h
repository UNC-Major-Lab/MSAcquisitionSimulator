//
// Created by Dennis Goldfarb on 6/27/15.
//

#ifndef MSACQUISITIONSIMULATOR_PROTEIN_H
#define MSACQUISITIONSIMULATOR_PROTEIN_H


#include <iostream>
#include <map>

class Protein {
public :

	Protein() {}

	Protein(std::string name, std::string sequence, std::string description, double abundance)
			: name(name), sequence(sequence), description(description), abundance(abundance) {//, residue2PTM() {

	}

	Protein(std::string name, std::string sequence, std::string description)
			: Protein(name, sequence, description, 0) {

	}

	~Protein() {}

	std::string name;
	std::string sequence;
	std::string description;
	double abundance;
};



#endif //MSACQUISITIONSIMULATOR_PROTEIN_H
