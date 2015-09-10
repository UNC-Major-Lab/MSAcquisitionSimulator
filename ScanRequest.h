//
// Created by Dennis Goldfarb on 8/20/15.
//

#ifndef MSACQUISITIONSIMULATOR_SCANREQUEST_H
#define MSACQUISITIONSIMULATOR_SCANREQUEST_H

#include "Peak.h"

class ScanRequest {
private:

public:
	ScanRequest(double min_mz, double max_mz, bool do_fragmentation, double max_injection_time, double target_total_ion_count) :
			min_mz(min_mz), max_mz(max_mz), do_fragmentation(do_fragmentation), max_injection_time(max_injection_time), target_total_ion_count(target_total_ion_count) {}

	double min_mz;
	double max_mz;
	double max_injection_time;
	double target_total_ion_count;
	bool do_fragmentation;
};


#endif //MSACQUISITIONSIMULATOR_SCANREQUEST_H
