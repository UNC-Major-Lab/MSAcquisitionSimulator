//
// Created by Dennis Goldfarb on 8/29/15.
//

#include "StaticExclusionListEntry.h"

bool StaticExclusionListEntry::is_match(StaticExclusionMatchQuery &q) const {
	return std::abs(q.mz - mz) <= q.mz_tolerance
		   && q.time >= start_time
		   && q.time <= end_time;
}

