//
// Created by Dennis Goldfarb on 8/20/15.
//

#ifndef MSACQUISITIONSIMULATOR_INSTRUMENT_H
#define MSACQUISITIONSIMULATOR_INSTRUMENT_H


#include <vector>
#include <algorithm>
#include "Peak.h"
#include "ScanRequest.h"

class Instrument {
private:

public:
	Instrument(double resolution, double dynamic_range, double ms1_scan_time, double ms2_scan_time, double scan_overhead_time,
			   double max_ms1_injection_time, double max_ms2_injection_time, double ms1_target_total_ion_count, double ms2_target_total_ion_count) :
			resolution(resolution), dynamic_range(dynamic_range), ms1_scan_time(ms1_scan_time), ms2_scan_time(ms2_scan_time),
			scan_overhead_time(scan_overhead_time), max_ms1_injection_time(max_ms1_injection_time), max_ms2_injection_time(max_ms2_injection_time),
			ms1_target_total_ion_count(ms1_target_total_ion_count), ms2_target_total_ion_count(ms2_target_total_ion_count) {}

	double resolution;
	double dynamic_range;
	double ms1_scan_time;
	double ms2_scan_time;
	double scan_overhead_time;
	double max_ms1_injection_time;
	double max_ms2_injection_time;
	double ms1_target_total_ion_count;
	double ms2_target_total_ion_count;

	void apply_dynamic_range(std::vector<BasicPeak> &peaks);
	double get_scan_overhead_time();
	double get_scan_acquisition_time(ScanRequest* scanRequest);
};


#endif //MSACQUISITIONSIMULATOR_INSTRUMENT_H
