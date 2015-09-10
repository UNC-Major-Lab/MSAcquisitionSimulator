//
// Created by Dennis Goldfarb on 8/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_RANDOMN_H
#define MSACQUISITIONSIMULATOR_RANDOMN_H


#include <vector>
#include <algorithm>
#include "Peak.h"
#include "AbstractTopN.h"
#include "TopNParameters.h"

class RandomN : public AbstractTopN {

private:
	std::list<BasicPeak> rank_peaks(std::vector<BasicPeak> &peaks);

public:

	RandomN(TopNParameters parameters) : AbstractTopN(parameters) {}

};


#endif //MSACQUISITIONSIMULATOR_RANDOMN_H
