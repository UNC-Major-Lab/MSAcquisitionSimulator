//
// Created by Dennis Goldfarb on 8/20/15.
//

#include "Instrument.h"

void Instrument::apply_dynamic_range(std::vector<BasicPeak> &peaks) {
	if (peaks.size() > 0) {
		double max_intensity = std::min_element(peaks.begin(), peaks.end(), BasicPeak::greater_intensity)->intensity;

		max_intensity /= dynamic_range;

		for (auto itr = peaks.begin(); itr != peaks.end();) {
			if (itr->intensity < max_intensity) {
				itr = peaks.erase(itr);
			} else {
				++itr;
			}
		}
	}
}

double Instrument::get_scan_acquisition_time(ScanRequest* scanRequest) {
	if (scanRequest->do_fragmentation) return ms2_scan_time;
	else return ms1_scan_time;
}

double Instrument::get_scan_overhead_time() {
	return scan_overhead_time;
}
