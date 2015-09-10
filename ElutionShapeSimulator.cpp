//
// Created by Dennis Goldfarb on 8/26/15.
//

#include "ElutionShapeSimulator.h"


double ElutionShapeSimulator::get_min_rt(double rt_center, double abundance) {

	double ln_ratio = log(abundance)-log(RT_PRUNE_THRESHOLD);
	double tau_ln_ratio = tau*ln_ratio;

	double min_rt = 0.5*(2*rt_center + tau_ln_ratio)
					- 0.5*(sqrt(8*sigma_squared*ln_ratio
						   + tau_ln_ratio*tau_ln_ratio));

	return min_rt;
}

double ElutionShapeSimulator::get_max_rt(double rt_center, double abundance) {
	double ln_ratio = log(abundance)-log(RT_PRUNE_THRESHOLD);
	double tau_ln_ratio = tau*ln_ratio;

	double max_rt = 0.5*(2*rt_center + tau_ln_ratio)
					+ 0.5*(sqrt(8*sigma_squared*ln_ratio
						   + tau_ln_ratio*tau_ln_ratio));

	return max_rt;
}

double ElutionShapeSimulator::get_abundance_over_time(double rt_center, double abundance, double start, double end) {
	// Numerical integration using Simpson's rule
	double area = abundance * ((end-start)/6) * (evaluate_EGH(rt_center, start)
									 + evaluate_EGH(rt_center, end)
									 + 4*evaluate_EGH(rt_center, (start+end)/2));

	return area;
}

double ElutionShapeSimulator::evaluate_EGH(double rt_center, double t) {
	double diff = t-rt_center;
	double diff_squared = diff*diff;
	double denom = 2*sigma_squared + tau*diff;

	if (denom > 0) {
		double x = exp(-diff_squared / denom);
		return x;
	}
	return 0;
}

void ElutionShapeSimulator::check_valid_parameters() {
	assert(sigma_squared>0);
	assert(tau>=0);
}

double ElutionShapeSimulator::get_max_abundance() {
	double AUC = 0;
	double diff = 0;
	double rt_start = 0;
	double step = .01;
	do {
		diff = get_abundance_over_time(0, 1, rt_start-step, rt_start);
		AUC += diff;
		rt_start-=step;

	} while (diff/AUC > 1e-10);

	rt_start = 0;
	do {
		diff = get_abundance_over_time(0, 1, rt_start, rt_start+step);
		AUC += diff;
		rt_start+=step;

	} while (diff/AUC > 1e-10);

	return AUC;
}
