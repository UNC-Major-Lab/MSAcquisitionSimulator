//
// Created by Dennis Goldfarb on 8/20/15.
//

#ifndef MSACQUISITIONSIMULATOR_SCAB_H
#define MSACQUISITIONSIMULATOR_SCAB_H


#include <vector>
#include "Peak.h"

class Scan {
private:

public:
	enum ScanType {MS1, MS2};

	Scan(std::vector<BasicPeak> peaks, double retention_time, double elapsed_time, int scan_id, ScanType scan_type) :
			peaks(peaks), retention_time(retention_time), elapsed_time(elapsed_time), scan_id(scan_id), scan_type(scan_type) {}

	std::vector<BasicPeak> peaks;
	double retention_time;
	double elapsed_time;
	int scan_id;
	ScanType scan_type;

};


#endif //MSACQUISITIONSIMULATOR_SCAB_H
