//
// Created by Dennis Goldfarb on 8/28/15.
//

#ifndef MSACQUISITIONSIMULATOR_STATICEXCLUSIONLIST_H
#define MSACQUISITIONSIMULATOR_STATICEXCLUSIONLIST_H

#include <vector>
#include "ExclusionList.h"
#include "StaticExclusionListEntry.h"
#include "Peak.h"

class StaticExclusionList : public ExclusionList<StaticExclusionListEntry> {

private:

public:
	StaticExclusionList(double mz_tolerance) : ExclusionList(mz_tolerance) {}

	void filter_by_list(std::vector<BasicPeak> &peaks, double current_time);
};


#endif //MSACQUISITIONSIMULATOR_STATICEXCLUSIONLIST_H
