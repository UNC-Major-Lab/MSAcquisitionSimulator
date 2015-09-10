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
	Instrument(double resolution, double dynamic_range) :
			resolution(resolution), dynamic_range(dynamic_range) {}

	double resolution;
	double dynamic_range;

	void apply_dynamic_range(std::vector<BasicPeak> &peaks);
	double get_scan_overhead_time();
	double get_scan_acquisition_time(ScanRequest* scanRequest);
};


#endif //MSACQUISITIONSIMULATOR_INSTRUMENT_H
