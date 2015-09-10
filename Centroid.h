//
// Created by Dennis Goldfarb on 8/30/15.
//

#ifndef MSACQUISITIONSIMULATOR_CENTROID_H
#define MSACQUISITIONSIMULATOR_CENTROID_H


class Centroid {
private:

public:
	Centroid(double mz, double mass, int charge, double intensity, int num_neutrons, double rt_center) :
			mz(mz), mass(mass), charge(charge), intensity(intensity), num_neutrons(num_neutrons), rt_center(rt_center) {}

	double mz;
	double mass;
	double rt_center;
	double intensity;
	int charge;
	int num_neutrons;

};


#endif //MSACQUISITIONSIMULATOR_CENTROID_H
