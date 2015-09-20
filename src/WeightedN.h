//
// Created by Dennis Goldfarb on 8/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_STOCHASTICN_H
#define MSACQUISITIONSIMULATOR_STOCHASTICN_H


#include <vector>
#include "Peak.h"
#include "AbstractTopN.h"
#include "TopNParameters.h"
#include "WeightedShuffle.h"

class WeightedN : public AbstractTopN {

private:
	std::list<BasicPeak> rank_peaks(std::vector<BasicPeak> &peaks);

public:
	WeightedN(const std::vector<std::string>& values) : AbstractTopN(values) {};

	WeightedN(TopNParameters parameters) : AbstractTopN(parameters) {}

};


#endif //MSACQUISITIONSIMULATOR_STOCHASTICN_H
