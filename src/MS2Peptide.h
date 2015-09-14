//
// Created by Dennis Goldfarb on 9/7/15.
//

#ifndef MSACQUISITIONSIMULATOR_MS2PEPTIDE_H
#define MSACQUISITIONSIMULATOR_MS2PEPTIDE_H


#include <string>

class MS2Peptide {
private:

public:
	MS2Peptide() {};

	MS2Peptide(double total_intensity, double precursor_ion_fraction, std::string modified_sequence) :
			total_intensity(total_intensity), precursor_ion_fraction(precursor_ion_fraction), modified_sequence(modified_sequence) {}

	double total_intensity;
	double precursor_ion_fraction;
	std::string modified_sequence;
};


#endif //MSACQUISITIONSIMULATOR_MS2PEPTIDE_H
