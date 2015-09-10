//
// Created by Dennis Goldfarb on 8/28/15.
//

#ifndef MSACQUISITIONSIMULATOR_TOPN_H
#define MSACQUISITIONSIMULATOR_TOPN_H

#include <iostream>
#include <vector>
#include <queue>
#include <list>
#include <algorithm>
#include "AcquisitionController.h"
#include "MS1Scan.h"
#include "MS2Scan.h"
#include "Scan.h"
#include "Peak.h"
#include "DynamicExclusionList.h"
#include "StaticExclusionList.h"
#include "SpectrumAnalyzer.h"
#include "MS2ScanRequest.h"
#include "TopNParameters.h"
#include "AbstractTopN.h"

class TopN : public AbstractTopN {

private:
	std::list<BasicPeak> rank_peaks(std::vector<BasicPeak> &peaks);

public:

	TopN(TopNParameters parameters) : AbstractTopN(parameters) {}

};


#endif //MSACQUISITIONSIMULATOR_TOPN_H
