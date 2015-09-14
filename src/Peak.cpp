//
// Created by Dennis Goldfarb on 8/28/15.
//

#include "Peak.h"

bool Peak::less_mz(const Peak &a, const Peak &b) {
	return a.mz < b.mz;
}

bool Peak::greater_intensity(const Peak &a, const Peak &b) {
	return a.intensity > b.intensity;
}
