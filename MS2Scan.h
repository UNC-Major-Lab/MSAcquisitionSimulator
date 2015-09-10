//
// Created by Dennis Goldfarb on 8/28/15.
//

#ifndef MSACQUISITIONSIMULATOR_MS2SCAN_H
#define MSACQUISITIONSIMULATOR_MS2SCAN_H

#include <string>
#include <map>
#include "Scan.h"
#include "Protein.h"

class MS2Scan : public Scan {
private:

public:
	MS2Scan(std::vector<BasicPeak> peaks, double retention_time, double elapsed_time, int scan_id,
			const ScanType &scan_type, BasicPeak precursor_peak, int parent_scan_id,
			std::map<std::string, double> peptide2intensity, double target_total_ion_count) :
			Scan(peaks, retention_time, elapsed_time, scan_id, scan_type), precursor_peak(precursor_peak),
			parent_scan_id(parent_scan_id), peptide2intensity(peptide2intensity), probability(0),
			target_total_ion_count(target_total_ion_count) {
		TIC = 0;
		for (auto itr = peptide2intensity.begin(); itr != peptide2intensity.end(); ++itr) {
			TIC+= itr->second;
		}
	}

	std::map<std::string, double> peptide2intensity;
	double TIC;
	BasicPeak precursor_peak;
	int parent_scan_id;
	std::string peptide;
	double probability;
	double target_total_ion_count;
	std::vector<Protein*> proteins;
};


#endif //MSACQUISITIONSIMULATOR_MS2SCAN_H
