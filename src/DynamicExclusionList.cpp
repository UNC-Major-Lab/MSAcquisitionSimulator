//
// Created by Dennis Goldfarb on 8/28/15.
//

#include "DynamicExclusionList.h"

void DynamicExclusionList::add_entry(double mz, double inserted_time) {
	entries.insert(DynamicExclusionListEntry(mz, inserted_time));
}

void DynamicExclusionList::filter_by_list(std::vector<BasicPeak> &peaks) {

	auto itr_entry = entries.begin();
	auto itr_peak = peaks.begin();

	for (; itr_entry != entries.end() && itr_peak != peaks.end();) {
		DynamicExclusionMatchQuery q = DynamicExclusionMatchQuery(itr_peak->mz, mz_tolerance);
		if (itr_entry->is_match(q)) {
			itr_peak = peaks.erase(itr_peak);
		} else if (itr_entry->mz < itr_peak->mz) {
			++itr_entry;
		} else {
			++itr_peak;
		}
	}
}

void DynamicExclusionList::remove_expired(double current_time) {
	for (auto itr = entries.begin(); itr != entries.end();) {
		if (current_time - itr->inserted_time >= exclusion_time) {
			itr = entries.erase(itr);
		} else {
			++itr;
		}
	}
}

void DynamicExclusionList::filter_clashing(std::list<BasicPeak> &peaks) {
	for (auto itr1 = peaks.begin(); itr1 != peaks.end(); ++itr1) {
		for (auto itr2 = std::next(itr1); itr2 != peaks.end();) {
			if (std::abs(itr1->mz - itr2->mz) <= mz_tolerance) {
				itr2 = peaks.erase(itr2);
			} else {
				++itr2;
			}
		}
	}
}
