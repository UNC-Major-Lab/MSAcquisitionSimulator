//
// Created by Dennis Goldfarb on 8/20/15.
//

#ifndef MSACQUISITIONSIMULATOR_ACQUISITIONCONTROLLER_H
#define MSACQUISITIONSIMULATOR_ACQUISITIONCONTROLLER_H

#include <memory>
#include "Scan.h"
#include "ScanRequest.h"
#include "MS1Scan.h"
#include "MS2Scan.h"

class AcquisitionController {
private:

public:

	virtual void process_scan(Scan* scan);
	//virtual void process_scan(MS1Scan scan);
	//virtual void process_scan(MS2Scan scan);

	virtual std::unique_ptr<ScanRequest> get_scan_request(double current_time);
};


#endif //MSACQUISITIONSIMULATOR_ACQUISITIONCONTROLLER_H
