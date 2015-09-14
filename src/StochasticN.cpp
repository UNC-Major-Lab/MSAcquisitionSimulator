//
// Created by Dennis Goldfarb on 8/29/15.
//

#include "StochasticN.h"

std::list<BasicPeak> StochasticN::rank_peaks(std::vector<BasicPeak> &peaks) {
	std::vector<double> weights;
	for (BasicPeak &p : peaks) weights.push_back(p.intensity);
	std::vector<BasicPeak> shuffled_peaks = WeightedShuffle::shuffle<BasicPeak>(peaks, weights);

	std::list<BasicPeak> targets;
	for (BasicPeak &p : shuffled_peaks) targets.push_back(p);

	return targets;
}
