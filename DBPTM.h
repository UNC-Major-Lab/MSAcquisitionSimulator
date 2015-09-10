//
// Created by Dennis Goldfarb on 9/6/15.
//

#ifndef MSACQUISITIONSIMULATOR_DBPTM_H
#define MSACQUISITIONSIMULATOR_DBPTM_H


#include "PTM.h"

class DBPTM : public PTM {
private:

public:
	DBPTM() {}

	DBPTM(const std::string name, const std::string abbreviation, const std::string residues, const std::string formula, bool is_static) :
			PTM(name, abbreviation, residues, formula), is_static(is_static) {
	}

	DBPTM& operator=(const DBPTM& other) {
		return *this;
	};

	bool is_static;
};


#endif //MSACQUISITIONSIMULATOR_DBPTM_H
