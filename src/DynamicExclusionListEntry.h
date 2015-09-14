//
// Created by Dennis Goldfarb on 8/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_DYNAMICEXCLUSIONLISTENTRY_H
#define MSACQUISITIONSIMULATOR_DYNAMICEXCLUSIONLISTENTRY_H

#include <cmath>
#include "DynamicExclusionMatchQuery.h"

class DynamicExclusionListEntry {
private:

public:

	DynamicExclusionListEntry(double mz, double inserted_time) : mz(mz), inserted_time(inserted_time) {}

	double mz;
	double inserted_time;

	bool is_match(DynamicExclusionMatchQuery & q) const;
	friend bool operator< (const DynamicExclusionListEntry& a, const DynamicExclusionListEntry& b) {
		return a.mz < b.mz;
	}
};


#endif //MSACQUISITIONSIMULATOR_DYNAMICEXCLUSIONLISTENTRY_H
