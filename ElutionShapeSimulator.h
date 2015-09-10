//
// Created by Dennis Goldfarb on 8/26/15.
//

#ifndef MSACQUISITIONSIMULATOR_ELUTIONSHAPESIMULATOR_H
#define MSACQUISITIONSIMULATOR_ELUTIONSHAPESIMULATOR_H

#include <cmath>
#include "Globals.h"

class ElutionShapeSimulator {
private:
	double evaluate_EGH(double rt_center, double t);
	void check_valid_parameters();

public:
	ElutionShapeSimulator(double tau, double sigma) : tau(tau), sigma_squared(sigma*sigma) {
		check_valid_parameters();
		normalization_factor = get_max_abundance();
	};

	double tau;
	double sigma_squared;
	double RT_PRUNE_THRESHOLD = std::max(1.0,PRUNE_THRESHOLD/2);
	double normalization_factor;

	double get_min_rt(double rt_center, double abundance);
	double get_max_rt(double rt_center, double abundance);
	double get_abundance_over_time(double rt_center, double abundance, double start, double end);
	double get_max_abundance();
};


#endif //MSACQUISITIONSIMULATOR_ELUTIONSHAPESIMULATOR_H
