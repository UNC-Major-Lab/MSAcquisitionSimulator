//
// Created by Dennis Goldfarb on 8/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_MS2SCANREQUEST_H
#define MSACQUISITIONSIMULATOR_MS2SCANREQUEST_H


#include "ScanRequest.h"
#include "Peak.h"

class MS2ScanRequest : public ScanRequest {
private:

public:
	MS2ScanRequest(BasicPeak peak, int parent_scan_id, double min_mz, double max_mz, double max_ms2_fill_time, double target_total_ion_count) :
			ScanRequest(min_mz, max_mz, true, max_ms2_fill_time, target_total_ion_count), peak(peak), parent_scan_id(parent_scan_id) {}

	BasicPeak peak;
	int parent_scan_id;
};


#endif //MSACQUISITIONSIMULATOR_MS2SCANREQUEST_H
