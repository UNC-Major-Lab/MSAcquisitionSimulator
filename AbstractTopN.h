//
// Created by Dennis Goldfarb on 8/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_ABSTRACTTOPN_H
#define MSACQUISITIONSIMULATOR_ABSTRACTTOPN_H


#include <vector>
#include <queue>
#include <list>
#include "ScanRequest.h"
#include "Peak.h"
#include "TopNParameters.h"
#include "SpectrumAnalyzer.h"
#include "DynamicExclusionList.h"
#include "StaticExclusionList.h"
#include "AcquisitionController.h"
#include "MS1Scan.h"
#include "MS2Scan.h"
#include "MS2ScanRequest.h"

class AbstractTopN : public AcquisitionController {

private:
	int last_ms1_scan_id;
	int num_ms2_since_ms1 = 0;
	std::list<BasicPeak> targets;
	StaticExclusionList static_exclusion_list;
	DynamicExclusionList dynamic_exclusion_list;
	SpectrumAnalyzer sa;
	TopNParameters parameters;

	virtual std::list<BasicPeak> rank_peaks(std::vector<BasicPeak> &peaks) = 0;

public:

	AbstractTopN(TopNParameters parameters) : parameters(parameters), static_exclusion_list(parameters.exclusion_list_mz_tolerance),
	dynamic_exclusion_list(parameters.exclusion_list_mz_tolerance, parameters.exclusion_time) {}

	void process_scan(Scan* scan);
	void process_scan(MS1Scan scan);
	void process_scan(MS2Scan scan);


	std::unique_ptr<ScanRequest> get_scan_request(double current_time);


};


#endif //MSACQUISITIONSIMULATOR_ABSTRACTTOPN_H
