//
// Created by Dennis Goldfarb on 8/20/15.
//

#include "AcquisitionController.h"

void AcquisitionController::process_scan(Scan* scan) {

}

std::unique_ptr<ScanRequest> AcquisitionController::get_scan_request(double current_time) {
	return std::unique_ptr<ScanRequest>(new ScanRequest(0, 0, false, 0, 0));
}
