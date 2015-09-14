//
// Created by Dennis Goldfarb on 8/20/15.
//


#include "Oracle.h"
#include "SpectrumAnalyzer.h"

std::unique_ptr<Scan> Oracle::get_scan_data(ScanRequest *scan_request, double time) {
	current_scan_id++;

	if (scan_request->do_fragmentation) {
		return get_ms2_data(static_cast<MS2ScanRequest*>(scan_request), time);
	} else {
		return get_ms1_data(scan_request, time);
	}
}

std::unique_ptr<MS1Scan> Oracle::get_ms1_data(ScanRequest* scan_request, double time) {
	std::vector<IonDAO*> ions = db->get_ions_at_rt(scan_request->min_mz, scan_request->max_mz, time);
	double elapsed_time = 0;

	std::vector<Peak> peaks = generate_peak_list(ions, time, elapsed_time, instrument->max_ms1_injection_time, instrument->ms1_target_total_ion_count);

	std::vector<BasicPeak> ms1_signals = generate_profile_MS_signals(peaks, scan_request->min_mz, scan_request->max_mz, instrument->resolution);
	ms1_signals = centroid_MS_signals(ms1_signals);
	instrument->apply_dynamic_range(ms1_signals);

	elapsed_time = std::max(elapsed_time, instrument->get_scan_acquisition_time(scan_request));
	elapsed_time += instrument->get_scan_overhead_time();

	std::sort(ms1_signals.begin(), ms1_signals.end(), BasicPeak::less_mz);
	return std::unique_ptr<MS1Scan>(new MS1Scan(ms1_signals, time, elapsed_time, current_scan_id, Scan::ScanType::MS1));
}



std::unique_ptr<MS2Scan> Oracle::get_ms2_data(MS2ScanRequest* scan_request, double time) {
	std::vector<IonDAO*> ions = db->get_ions_at_rt(scan_request->min_mz, scan_request->max_mz, time);
	double elapsed_time = 0;

	std::map<std::string, double> peptide2intensity = calculate_peptide_abundance(ions, time, elapsed_time, instrument->max_ms2_injection_time, instrument->ms2_target_total_ion_count);
	std::vector<Peak> peaks = generate_peak_list(ions, time, elapsed_time, instrument->max_ms2_injection_time, instrument->ms2_target_total_ion_count);

	std::vector<BasicPeak> ms2_signals = generate_profile_MS_signals(peaks, scan_request->min_mz, scan_request->max_mz, instrument->resolution);
	ms2_signals = centroid_MS_signals(ms2_signals);

	elapsed_time = std::max(elapsed_time, instrument->get_scan_acquisition_time(scan_request));
	elapsed_time += instrument->get_scan_overhead_time();

	return std::unique_ptr<MS2Scan>(new MS2Scan(ms2_signals, time, elapsed_time, current_scan_id, Scan::ScanType::MS2, scan_request->peak, scan_request->parent_scan_id, peptide2intensity, instrument->ms2_target_total_ion_count));
}



std::map<std::string, double> Oracle::calculate_peptide_abundance(std::vector<IonDAO *> &ions, double current_time,
																  double &elapsed_time, double max_injection_time,
																  double target_total_ion_count) {
	std::map<std::string, double> peptide2intensity;

	double TIC = 0;
	double step_size = 0.001;

	for (IonDAO *ion : ions) {
		if (peptide2intensity.find(ion->modified_sequence) == peptide2intensity.end()) {
			peptide2intensity[ion->modified_sequence] = 0;
		}
	}

	for (elapsed_time = 0; elapsed_time < max_injection_time && TIC < target_total_ion_count; elapsed_time+=step_size) {
		for (IonDAO *ion : ions) {
			double abundance = elution_shape_simulator.get_abundance_over_time(ion->rt, ion->abundance, current_time+elapsed_time, current_time+elapsed_time+step_size);
			peptide2intensity[ion->modified_sequence]+= abundance;
			TIC += abundance;
		}
	}
	for (auto peptide = peptide2intensity.begin(); peptide != peptide2intensity.end();) {
		if (peptide->second <= 0) {
			peptide = peptide2intensity.erase(peptide);
		} else {
			++peptide;
		}
	}

	return peptide2intensity;
}

std::vector<Peak> Oracle::generate_peak_list(std::vector<IonDAO*> &ions, double current_time, double &elapsed_time, double max_injection_time, double target_total_ion_count) {
	// integrate abundance until target ions reached

	double TIC = 0;
	double step_size = 0.005;
	std::vector<Peak> peaks;
	for (IonDAO *ion : ions) peaks.push_back(Peak(ion->mz, ion->mass, ion->charge, 0, ion->neutrons));

	for (elapsed_time = 0; elapsed_time < max_injection_time && TIC < target_total_ion_count; elapsed_time+=step_size) {
		for (int i = 0; i < ions.size(); i++) {
			IonDAO *ion = ions[i];
			double abundance = elution_shape_simulator.get_abundance_over_time(ion->rt, ion->abundance, current_time+elapsed_time, current_time+elapsed_time+step_size);
			peaks[i].intensity += abundance;
			TIC += abundance;
		}
	}
	for (auto peak = peaks.begin(); peak != peaks.end();) {
		if (peak->intensity <= 0) {
			peak = peaks.erase(peak);
		} else {
			++peak;
		}
	}

	return peaks;
}

std::vector<BasicPeak> Oracle::generate_profile_MS_signals(std::vector<Peak> &peaks, double min_mz, double max_mz, double resolution) {
	std::vector<BasicPeak> ms_signals;

	double max_mz_distance = 1000/resolution;


	std::sort(peaks.begin(), peaks.end(), BasicPeak::less_mz);
	double step_size = 0.01;
	int step = 0;
	int min_peak = 0;

	for (double current_mz = min_mz; current_mz <= max_mz; current_mz += step_size) {
		ms_signals.push_back(BasicPeak(current_mz, 0));

		for (int i = min_peak; i < peaks.size(); i++) {
			Peak &peak = peaks[i];
			if (current_mz - peak.mz > max_mz_distance) {
				min_peak = i;
				continue;
			} else if (peak.mz - current_mz > max_mz_distance) {
				break;
			} else {
				double fwhm = peak.mz / resolution;
				double signal_intensity = peak.intensity * pdf_lorentzian(fwhm, current_mz, peak.mz);
				//std::cout << fwhm << "\t" << peak.intensity << "\t" << pdf_lorentzian(fwhm, current_mz, peak.mz) << "\t" << signal_intensity << "\t" << current_mz << "\t" << peak.mz << "\t" << std::endl;
				ms_signals[step].intensity += signal_intensity;
			}
		}
		++step;
	}

	return ms_signals;
}

std::vector<BasicPeak> Oracle::centroid_MS_signals(std::vector<BasicPeak> &profile_signals) {
	std::vector<BasicPeak> centroided_signals;


	// special case for left side
	if (profile_signals.size() > 1 && profile_signals[0].intensity > profile_signals[1].intensity) {
		centroided_signals.push_back(profile_signals[0]);
	}

	for (int i = 1; i < profile_signals.size()-1; i++) {
		double intensity = profile_signals[i].intensity;
		double left_intensity = profile_signals[i-1].intensity;
		double right_intensity = profile_signals[i+1].intensity;

		if (intensity > left_intensity && intensity > right_intensity) {
			centroided_signals.push_back(profile_signals[i]);
		}
	}

	// special case for right side
	if (profile_signals.size() > 1 && profile_signals[profile_signals.size()-1].intensity > profile_signals[profile_signals.size()-2].intensity) {
		centroided_signals.push_back(profile_signals[profile_signals.size()-1]);
	}

	return centroided_signals;
}

double Oracle::pdf_lorentzian(double fwhm, double x, double x0) {
	boost::math::cauchy_distribution<double> cauchy(0., fwhm / 2.0);
	double pdf = (boost::math::pdf(cauchy, (x-.005)-x0) + boost::math::pdf(cauchy, (x+.005)-x0))/2;
	return pdf*.01;
}


