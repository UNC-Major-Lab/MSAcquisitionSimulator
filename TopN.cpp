//
// Created by Dennis Goldfarb on 8/28/15.
//

#include "TopN.h"

std::list<BasicPeak> TopN::rank_peaks(std::vector<BasicPeak> &peaks) {
	std::list<BasicPeak> targets;

	sort(peaks.begin(), peaks.end(), BasicPeak::greater_intensity);

	for (BasicPeak &p : peaks) targets.push_back(p);

	return targets;
}
