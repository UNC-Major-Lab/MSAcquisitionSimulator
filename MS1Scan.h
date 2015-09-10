//
// Created by Dennis Goldfarb on 8/28/15.
//

#ifndef MSACQUISITIONSIMULATOR_MS1SCAN_H
#define MSACQUISITIONSIMULATOR_MS1SCAN_H


#include "Scan.h"

class MS1Scan : public Scan {

public:
	MS1Scan(std::vector<BasicPeak> peaks, double retention_time, double elapsed_time, int scan_id,
			const ScanType &scan_type) : Scan(peaks, retention_time, elapsed_time, scan_id, scan_type) { }
};


#endif //MSACQUISITIONSIMULATOR_MS1SCAN_H
