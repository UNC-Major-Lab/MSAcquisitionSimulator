//
// Created by Dennis Goldfarb on 8/28/15.
//

#include "StaticExclusionList.h"

void StaticExclusionList::filter_by_list(std::vector<BasicPeak> &peaks,double current_time) {

	auto itr_entry = entries.begin();
	auto itr_peak = peaks.begin();

	for (; itr_entry != entries.end() && itr_peak != peaks.end();) {
		StaticExclusionMatchQuery q = StaticExclusionMatchQuery(itr_peak->mz, mz_tolerance, current_time);
		if (itr_entry->is_match(q)) {
			itr_peak = peaks.erase(itr_peak);
		} else if (itr_entry->mz < itr_peak->mz) {
			++itr_entry;
		} else {
			++itr_peak;
		}
	}

}
