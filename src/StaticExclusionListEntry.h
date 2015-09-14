//
// Created by Dennis Goldfarb on 8/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_STATICEXCLUSIONLISTENTRY_H
#define MSACQUISITIONSIMULATOR_STATICEXCLUSIONLISTENTRY_H

#include <cmath>
#include "StaticExclusionMatchQuery.h"

class StaticExclusionListEntry {
private:

public:

	StaticExclusionListEntry(double mz, double start_time, double end_time) :
			mz(mz), start_time(start_time), end_time(end_time) {}

	bool is_match(StaticExclusionMatchQuery &q) const;
	friend bool operator< (const StaticExclusionMatchQuery& a, const StaticExclusionMatchQuery& b) {
		return a.mz < b.mz;
	}

	double mz;
	double start_time;
	double end_time;
};


#endif //MSACQUISITIONSIMULATOR_STATICEXCLUSIONLISTENTRY_H
