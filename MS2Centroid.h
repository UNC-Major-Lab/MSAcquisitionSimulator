//
// Created by Dennis Goldfarb on 9/7/15.
//

#ifndef MSACQUISITIONSIMULATOR_MS2CENTROID_H
#define MSACQUISITIONSIMULATOR_MS2CENTROID_H

#include <string>
#include "Centroid.h"


class MS2Centroid : public Centroid {
private:

public:
	MS2Centroid(double mz, double mass, int charge, double intensity, int num_neutrons, double rt_center, std::string modified_sequence) :
			Centroid(mz, mass, charge, intensity, num_neutrons, rt_center), modified_sequence(modified_sequence) {};

	std::string modified_sequence;
};


#endif //MSACQUISITIONSIMULATOR_MS2CENTROID_H
