//
// Created by Dennis Goldfarb on 8/20/15.
//

#ifndef MSACQUISITIONSIMULATOR_ORACLE_H
#define MSACQUISITIONSIMULATOR_ORACLE_H


#include <cmath>
#include <boost/math/distributions.hpp>
#include "ScanRequest.h"
#include "Scan.h"
#include "Instrument.h"
#include "MS1Scan.h"
#include "MS2Scan.h"
#include "MS2ScanRequest.h"
#include "GroundTruthText.h"

class Oracle {
private:
	Instrument* instrument;
	GroundTruthText* db;
	ElutionShapeSimulator elution_shape_simulator;
	int current_scan_id;

	std::unique_ptr<MS1Scan> get_ms1_data(ScanRequest* scan_request, double time);
	std::unique_ptr<MS2Scan> get_ms2_data(MS2ScanRequest* scan_request, double time);


	std::map<std::string, double> calculate_peptide_abundance(std::vector<IonDAO*> &ions, double current_time, double &elapsed_time, double max_injection_time, double target_total_ion_count);
	std::vector<Peak> generate_peak_list(std::vector<IonDAO*> &ions, double current_time, double &elapsed_time, double max_injection_time, double target_total_ion_count);

	std::vector<BasicPeak> generate_profile_MS_signals(std::vector<Peak> &peaks, double min_mz, double max_mz, double resolution);
	std::vector<BasicPeak> centroid_MS_signals(std::vector<BasicPeak> &profile_signals);
	double pdf_lorentzian(double fwhm, double x, double x0);


public:
	Oracle(GroundTruthText* db, Instrument* instrument, ElutionShapeSimulator elution_shape_simulator) :
			db(db), instrument(instrument), elution_shape_simulator(elution_shape_simulator) {
		current_scan_id = 0;
	};

	std::unique_ptr<Scan> get_scan_data(ScanRequest* scan_request, double time);
};


#endif //MSACQUISITIONSIMULATOR_ORACLE_H
