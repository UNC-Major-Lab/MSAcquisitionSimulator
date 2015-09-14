//
// Created by Dennis Goldfarb on 9/6/15.
//

#ifndef MSACQUISITIONSIMULATOR_DBENZYME_H
#define MSACQUISITIONSIMULATOR_DBENZYME_H

#include <string>
#include "Enzyme.h"
#include "Globals.h"

class DBEnzyme : public Enzyme {
private:

public:
	DBEnzyme() {}
	DBEnzyme(std::string name, std::string residues, std::string blocking_residues, TERMINUS terminus) :
			Enzyme(name), residues(residues), blocking_residues(blocking_residues), terminus(terminus) {}

	std::string residues;
	std::string blocking_residues;
	TERMINUS terminus;

	DBEnzyme& operator=(const DBEnzyme& other) {
		return *this;
	};
};


#endif //MSACQUISITIONSIMULATOR_DBENZYME_H
