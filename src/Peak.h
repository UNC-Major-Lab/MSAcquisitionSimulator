//
// Created by Dennis Goldfarb on 8/28/15.
//

#ifndef MSACQUISITIONSIMULATOR_PEAK_H
#define MSACQUISITIONSIMULATOR_PEAK_H


#include "BasicPeak.h"

class Peak : public BasicPeak {
private:

public:
	Peak(double mz, double mass, int charge, double intensity, int num_neutrons) :
			BasicPeak(mz, intensity), mass(mass), charge(charge), num_neutrons(num_neutrons) {}

	double mass;
	int charge;
	int num_neutrons;

	static bool greater_intensity(const Peak& a, const Peak& b);
	static bool less_mz(const Peak& a, const Peak& b);
};


#endif //MSACQUISITIONSIMULATOR_PEAK_H
