//
// Created by Dennis Goldfarb on 9/3/15.
//

#include "BasicPeak.h"

bool BasicPeak::less_mz(const BasicPeak &a, const BasicPeak &b) {
	return a.mz < b.mz;
}

bool BasicPeak::greater_intensity(const BasicPeak &a, const BasicPeak &b) {
	return a.intensity > b.intensity;
}