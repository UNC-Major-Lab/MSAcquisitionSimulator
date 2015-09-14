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
	if (parameters.dynamic_exclusion_enabled) {
		dynamic_exclusion_list.remove_expired(scan.retention_time);
	}

	sa.deisotope(scan.peaks);

	if (parameters.dynamic_exclusion_enabled) {
		static_exclusion_list.filter_by_list(scan.peaks, scan.retention_time);
		dynamic_exclusion_list.filter_by_list(scan.peaks);
	}

	targets = rank_peaks(scan.peaks);

	if (parameters.dynamic_exclusion_enabled) {
		dynamic_exclusion_list.filter_clashing(targets);
	}

	last_ms1_scan_id = scan.scan_id;
}

void AbstractTopN::process_scan(MS2Scan scan) {
	if (parameters.dynamic_exclusion_enabled) {
		dynamic_exclusion_list.add_entry(scan.precursor_peak.mz, scan.retention_time);
	}
}

std::unique_ptr<ScanRequest> AbstractTopN::get_scan_request(double current_time) {

	if (num_ms2_since_ms1 < parameters.n && targets.size() > 0) {
		num_ms2_since_ms1++;

		BasicPeak target = targets.front();
		targets.pop_front();

		double min_mz = target.mz - parameters.isolation_width;
		double max_mz = target.mz + parameters.isolation_width;
		return std::unique_ptr<MS2ScanRequest>(new MS2ScanRequest(target, last_ms1_scan_id, min_mz, max_mz));
	} else {
		num_ms2_since_ms1=0;
		return std::unique_ptr<ScanRequest>(new ScanRequest(parameters.ms1_min_mz, parameters.ms1_max_mz, false));
	}
}


void AbstractTopN::validate(const std::vector<std::string> &values) {
	int num_ms2;
	double ms1_min_mz;
	double ms1_max_mz;
	double ms2_isolation_width;
	double dynamic_exclusion_tolerance;
	int dynamic_exclusion_time;
	bool dynamic_exclusion_enabled;

	typedef boost::escaped_list_separator<char> separator_type;
	separator_type separator("\\",    // The escape characters.
							 "\t ",    // The separator characters.
							 "\"\'"); // The quote characters.

	boost::program_options::options_description general("AbstractTopN");
	general.add_options()
			("num_ms2", boost::program_options::value<int>(&num_ms2)->default_value(10), "num_ms2")
			("ms1_min_mz", boost::program_options::value<double>(&ms1_min_mz)->default_value(200), "ms1_min_mz")
			("ms1_max_mz", boost::program_options::value<double>(&ms1_max_mz)->default_value(3000), "ms1_max_mz")
			("ms2_isolation_width", boost::program_options::value<double>(&ms2_isolation_width)->default_value(1), "ms2_isolation_width")
			("dynamic_exclusion_enabled", boost::program_options::value<bool>(&dynamic_exclusion_enabled)->default_value(true), "dynamic_exclusion_enabled")
			("dynamic_exclusion_tolerance", boost::program_options::value<double>(&dynamic_exclusion_tolerance)->default_value(0.05), "dynamic_exclusion_tolerance")
			("dynamic_exclusion_time", boost::program_options::value<int>(&dynamic_exclusion_time)->default_value(30), "dynamic_exclusion_time")
			;

	boost::tokenizer<separator_type> tokens(values[0], separator);
	std::vector<std::string> result;

	std::copy(tokens.begin(), tokens.end(), std::back_inserter(result));

	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::command_line_parser(result).options(general).run(), vm);
	boost::program_options::notify(vm);

	parameters = TopNParameters(ms2_isolation_width, ms1_min_mz, ms1_max_mz, dynamic_exclusion_tolerance,
								dynamic_exclusion_time, num_ms2, dynamic_exclusion_enabled);
	static_exclusion_list = StaticExclusionList(dynamic_exclusion_tolerance);
	dynamic_exclusion_list = DynamicExclusionList(dynamic_exclusion_tolerance, dynamic_exclusion_time);
}
