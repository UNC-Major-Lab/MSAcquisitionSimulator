//
// Created by Dennis Goldfarb on 8/29/15.
//

#include "AbstractTopN.h"

void AbstractTopN::process_scan(Scan* scan) {
	switch (scan->scan_type) {
		case Scan::ScanType::MS1:
			process_scan(*static_cast<MS1Scan*>(scan));
			break;
		case Scan::ScanType::MS2:
			process_scan(*static_cast<MS2Scan*>(scan));
	}
}

void AbstractTopN::process_scan(MS1Scan scan) {
	dynamic_exclusion_list.remove_expired(scan.retention_time);

	sa.deisotope(scan.peaks);

	static_exclusion_list.filter_by_list(scan.peaks, scan.retention_time);
	dynamic_exclusion_list.filter_by_list(scan.peaks);

	targets = rank_peaks(scan.peaks);

	dynamic_exclusion_list.filter_clashing(targets);

	last_ms1_scan_id = scan.scan_id;
}

void AbstractTopN::process_scan(MS2Scan scan) {
	dynamic_exclusion_list.add_entry(scan.precursor_peak.mz, scan.retention_time);
}

std::unique_ptr<ScanRequest> AbstractTopN::get_scan_request(double current_time) {

	if (num_ms2_since_ms1 < parameters.n && targets.size() > 0) {
		num_ms2_since_ms1++;

		BasicPeak target = targets.front();
		targets.pop_front();

		double min_mz = target.mz - parameters.isolation_width;
		double max_mz = target.mz + parameters.isolation_width;
		return std::unique_ptr<MS2ScanRequest>(new MS2ScanRequest(target, last_ms1_scan_id, min_mz, max_mz, parameters.max_ms2_injection_time, parameters.ms2_target_total_ion_count));
	} else {
		num_ms2_since_ms1=0;
		return std::unique_ptr<ScanRequest>(new ScanRequest(parameters.ms1_min_mz, parameters.ms1_max_mz, false, parameters.max_ms1_injection_time, parameters.ms1_target_total_ion_count));
	}
}


