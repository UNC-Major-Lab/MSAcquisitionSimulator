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
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/bind.hpp>

/*
 * This is the inferface that all new acquisition algorithms must inherit
 */
class AcquisitionController {
private:

public:

	virtual void process_scan(Scan* scan);

	virtual std::unique_ptr<ScanRequest> get_scan_request(double current_time);

	virtual void validate(const std::vector<std::string>& values);
};


#endif //MSACQUISITIONSIMULATOR_ACQUISITIONCONTROLLER_H
