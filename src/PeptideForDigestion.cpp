//
// Created by Dennis Goldfarb on 7/13/15.
//

#include "PeptideForDigestion.h"

const char& PeptideForDigestion::operator[](std::size_t i) const {
	return protein->sequence[i];
}
