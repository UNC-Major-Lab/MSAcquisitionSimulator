//
// Created by Dennis Goldfarb on 8/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_TOPNPARAMETERS_H
#define MSACQUISITIONSIMULATOR_TOPNPARAMETERS_H


class TopNParameters {
private:

public:
	TopNParameters(double isolation_width, double ms1_min_mz, double ms1_max_mz, double exclusion_list_mz_tolerance,
				   double exclusion_time, int n, bool dynamic_exclusion_enabled, double max_ms1_injection_time, double max_ms2_injection_time,
				   double ms1_target_total_ion_count, double ms2_target_total_ion_count) :
			isolation_width(isolation_width), ms1_max_mz(ms1_max_mz),
			ms1_min_mz(ms1_min_mz), exclusion_list_mz_tolerance(exclusion_list_mz_tolerance),
			exclusion_time(exclusion_time), n(n), dynamic_exclusion_enabled(dynamic_exclusion_enabled),
			max_ms1_injection_time(max_ms1_injection_time), max_ms2_injection_time(max_ms2_injection_time),
			ms1_target_total_ion_count(ms1_target_total_ion_count), ms2_target_total_ion_count(ms2_target_total_ion_count) {}

	double isolation_width;
	double ms1_min_mz;
	double ms1_max_mz;
	double exclusion_list_mz_tolerance;
	double exclusion_time;
	double max_ms1_injection_time;
	double max_ms2_injection_time;
	double ms1_target_total_ion_count;
	double ms2_target_total_ion_count;
	int n;
	bool dynamic_exclusion_enabled;
};


#endif //MSACQUISITIONSIMULATOR_TOPNPARAMETERS_H
