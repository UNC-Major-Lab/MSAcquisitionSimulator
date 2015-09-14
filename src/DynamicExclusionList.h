//
// Created by Dennis Goldfarb on 8/28/15.
//

#ifndef MSACQUISITIONSIMULATOR_DYNAMICEXCLUSIONLIST_H
#define MSACQUISITIONSIMULATOR_DYNAMICEXCLUSIONLIST_H

#include <iostream>
#include <list>
#include <vector>
#include "Peak.h"
#include "DynamicExclusionListEntry.h"
#include "DynamicExclusionMatchQuery.h"
#include "ExclusionList.h"

class DynamicExclusionList : public ExclusionList<DynamicExclusionListEntry> {

private:

public:
	DynamicExclusionList() {}
	DynamicExclusionList(double mz_tolerance, double exclusion_time) : ExclusionList(mz_tolerance), exclusion_time(exclusion_time) {}

	double exclusion_time;

	void add_entry(double mz, double inserted_time);

	void filter_by_list(std::vector<BasicPeak> &peaks);
	void remove_expired(double current_time);
	void filter_clashing(std::list<BasicPeak> &peaks);
};


#endif //MSACQUISITIONSIMULATOR_DYNAMICEXCLUSIONLIST_H
