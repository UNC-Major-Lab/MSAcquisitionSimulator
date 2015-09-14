//
// Created by Dennis Goldfarb on 8/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_TOPNPARAMETERS_H
#define MSACQUISITIONSIMULATOR_TOPNPARAMETERS_H


class TopNParameters {
private:

public:
	TopNParameters() {}

	TopNParameters(double isolation_width, double ms1_min_mz, double ms1_max_mz, double exclusion_list_mz_tolerance,
				   double exclusion_time, int n, bool dynamic_exclusion_enabled) :
			isolation_width(isolation_width), ms1_max_mz(ms1_max_mz),
			ms1_min_mz(ms1_min_mz), exclusion_list_mz_tolerance(exclusion_list_mz_tolerance),
			exclusion_time(exclusion_time), n(n), dynamic_exclusion_enabled(dynamic_exclusion_enabled) {}

	double isolation_width;
	double ms1_min_mz;
	double ms1_max_mz;
	double exclusion_list_mz_tolerance;
	double exclusion_time;
	int n;
	bool dynamic_exclusion_enabled;
};


#endif //MSACQUISITIONSIMULATOR_TOPNPARAMETERS_H
