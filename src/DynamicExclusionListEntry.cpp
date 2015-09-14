//
// Created by Dennis Goldfarb on 8/29/15.
//

#include "DynamicExclusionListEntry.h"

bool DynamicExclusionListEntry::is_match(DynamicExclusionMatchQuery &q) const {
	return std::abs(q.mz - mz) <= q.mz_tolerance;
}
