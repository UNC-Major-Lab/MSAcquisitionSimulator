//
// Created by Dennis Goldfarb on 8/29/15.
//

#include "RandomN.h"

std::list<BasicPeak> RandomN::rank_peaks(std::vector<BasicPeak> &peaks) {
	std::list<BasicPeak> targets;

	random_shuffle(peaks.begin(), peaks.end());

	for (BasicPeak &p : peaks) targets.push_back(p);

	return targets;
}
