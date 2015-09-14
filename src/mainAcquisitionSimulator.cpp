//
// Created by Dennis Goldfarb on 8/19/15.
//

#include <iostream>
#include <fstream>
#include <chrono>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/bind.hpp>

#include "Globals.h"
#include "AcquisitionController.h"
#include "Oracle.h"
#include "Sequencer.h"
#include "TopN.h"
#include "MzMLWriter.h"
#include "DBEnzyme.h"
#include "DBPTM.h"
#include "FidoWriter.h"
#include "GroundTruthText.h"
#include "StochasticN.h"
#include "RandomN.h"

#include "MSAcquisitionSimulatorConfig.h"

/*
 * Begin: Command line and config file parsing for user defined classes
 */
void validate(boost::any& v, const std::vector<std::string>& values, TERMINUS* target_type, int) {
	if (values[0][0] == 'N') v = TERMINUS::N_TERM;
	else if (values[0][0] == 'C') v = TERMINUS::C_TERM;
}

void validate(boost::any& v, const std::vector<std::string>& values, DBEnzyme* target_type, int) {
	std::string name;
	std::string residues;
	std::string blocking_residues;
	TERMINUS terminus;

	typedef boost::escaped_list_separator<char> separator_type;
	separator_type separator("\\",    // The escape characters.
							 "\t ",    // The separator characters.
							 "\"\'"); // The quote characters.

	boost::program_options::options_description general("Enzymes");
	general.add_options()
			("name,n", boost::program_options::value<std::string>(&name), "")
			("residues,r", boost::program_options::value<std::string>(&residues), "")
			("blocking_residues,b", boost::program_options::value<std::string>(&blocking_residues), "")
			("terminus,t", boost::program_options::value<TERMINUS>(&terminus), "")
			;

	boost::tokenizer<separator_type> tokens(values[0], separator);
	std::vector<std::string> result;

	std::copy(tokens.begin(), tokens.end(), std::back_inserter(result));

	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::command_line_parser(result).options(general).run(), vm);
	boost::program_options::notify(vm);

	std::cout << "Enzyme registered: " << name << " " << residues << " " << blocking_residues  << std::endl;

	v = DBEnzyme(name, residues, blocking_residues, terminus);
}

void validate(boost::any& v, const std::vector<std::string>& values, DBPTM* target_type, int) {
	std::string name;
	std::string abbreviation;
	std::string residues;
	std::string formula;
	bool is_static;

	typedef boost::escaped_list_separator<char> separator_type;
	separator_type separator("\\",    // The escape characters.
							 "\t ",    // The separator characters.
							 "\"\'"); // The quote characters.

	boost::program_options::options_description general("PTMs");
	general.add_options()
			("name,n", boost::program_options::value<std::string>(&name), "")
			("abbreviation,b", boost::program_options::value<std::string>(&abbreviation), "")
			("residues,r", boost::program_options::value<std::string>(&residues), "")
			("formula,m", boost::program_options::value<std::string>(&formula), "" )
			("is_static,s", boost::program_options::value<bool>(&is_static), "" )
			;

	boost::tokenizer<separator_type> tokens(values[0], separator);
	std::vector<std::string> result;

	std::copy(tokens.begin(), tokens.end(), std::back_inserter(result));

	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::command_line_parser(result).options(general).run(), vm);
	boost::program_options::notify(vm);

	std::cout << "PTM registered: " << name << " " << residues << " " << formula << std::endl;

	v = DBPTM(name, abbreviation, residues, formula, is_static);
}
/*
 * End: Command line and config file parsing for user defined classes
 */


/*
 * Simple map to convert algorithm name to the actual AcquisitionController.
 * You must register your new controller here.
 */
AcquisitionController* get_controller(std::string name, std::vector<std::string> acquisition_param_values) {
	if (name == "TopN") {
		return new TopN(acquisition_param_values);
	} else if (name == "RandomN") {
		return new RandomN(acquisition_param_values);
	} else if (name == "StochasticN") {
		return new StochasticN(acquisition_param_values);
	}

	std::cout << "Acquisition algorithm name not recognized: " << name << ". Defaulting to TopN.";
	return new TopN(acquisition_param_values); // DEFAULT
}



int main(int argc, const char ** argv) {

	std::cout << "MSAcquisitionSimulator Version " << MSAcquisitionSimulator_VERSION_MAJOR << "." << MSAcquisitionSimulator_VERSION_MINOR << "." << MSAcquisitionSimulator_VERSION_PATCH << std::endl;
	std::cout << "AcquisitionSimulator Version " << AcquisitionSimulator_VERSION_MAJOR << "." << AcquisitionSimulator_VERSION_MINOR << "." << AcquisitionSimulator_VERSION_PATCH << std::endl;

	std::string sqlite_in_path;
	std::string mzml_out_path;
	std::string param_file_path;
	std::string fido_out_path;
	std::string fasta_in_path;
	std::string acquisition_algorithm_name;

	std::vector<std::string> acquisition_param_values;

	double ms1_scan_time;
	double ms2_scan_time;
	double scan_overhead_time;
	int acquisition_length;

	double elution_tau;
	double elution_sigma;

	double resolution;
	double dynamic_range;

	double db_search_min_mass;
	double db_search_max_mass;
	double db_search_mass_tolerance;
	int db_search_max_missed_cleavages;
	int db_search_max_dynamic_mods;
	int db_search_min_enzymatic_termini;
	double null_lambda;

	double max_ms1_injection_time;
	double max_ms2_injection_time;
	double ms1_target_total_ion_count;
	double ms2_target_total_ion_count;

	std::vector<DBPTM> PTMs;
	std::vector<DBEnzyme> enzymes;


	//region Command line specification
	boost::program_options::options_description general("USAGE: AcquisitionSimulator [options] ground_truth.tab\n\nOptions");
	general.add_options()
			("help", "Print usage and exit.")
			("conf,c", boost::program_options::value<std::string>(&param_file_path)->default_value("acquisition.conf"), "Input path to config file.")
			("mzml_out_path,o", boost::program_options::value<std::string>(&mzml_out_path)->default_value("sample.mzML"), "output path for mzML file.")
			("fido_out_path,f", boost::program_options::value<std::string>(&fido_out_path)->default_value("sample.fido"), "output path for fido file.")
			;

	boost::program_options::options_description hidden("");
	hidden.add_options()
			("sqlite_in_path", boost::program_options::value<std::string>(&sqlite_in_path), "input path for ground truth file made by GroundTruthSimulator.")
			;

	boost::program_options::options_description all("Allowed options");
	all.add(general).add(hidden);

	boost::program_options::positional_options_description p;
	p.add("sqlite_in_path", -1);

	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(all).positional(p).run(), vm);
	boost::program_options::notify(vm);
	//endregion

	//region Command line processing
	if (vm.count("help")) {
		std::cout << general << std::endl;
		return 0;
	}

	//region config file specification
	boost::program_options::options_description config("Configuration file options");
	config.add_options()
			("acquisition_algorithm", boost::program_options::value<std::string>(&acquisition_algorithm_name)->default_value("TopN"), "acquisition algorithm")
			("acquisition_algorithm_params", boost::program_options::value<std::vector<std::string>>(&acquisition_param_values)->multitoken(), "acquisition_algorithm_params")
			("ms1_scan_time", boost::program_options::value<double>(&ms1_scan_time)->default_value(0.256), "ms1_scan_time")
			("ms2_scan_time", boost::program_options::value<double>(&ms2_scan_time)->default_value(0.064), "ms2_scan_time")
			("scan_overhead_time", boost::program_options::value<double>(&scan_overhead_time)->default_value(0.015), "scan_overhead_time")
			("acquisition_length", boost::program_options::value<int>(&acquisition_length)->default_value(3600), "acquisition_length")
			("fasta", boost::program_options::value<std::string>(&fasta_in_path), "acquisition fasta_in_path")
			("elution_tau", boost::program_options::value<double>(&elution_tau)->default_value(4), "elution_tau")
			("elution_sigma", boost::program_options::value<double>(&elution_sigma)->default_value(6), "elution_sigma")
			("resolution", boost::program_options::value<double>(&resolution)->default_value(60000), "resolution")
			("dynamic_range", boost::program_options::value<double>(&dynamic_range)->default_value(5000), "dynamic_range")
			("db_search_min_mass", boost::program_options::value<double>(&db_search_min_mass)->default_value(200), "db_search_min_mass")
			("db_search_max_mass", boost::program_options::value<double>(&db_search_max_mass)->default_value(9000), "db_search_max_mass")
			("db_search_max_missed_cleavages", boost::program_options::value<int>(&db_search_max_missed_cleavages)->default_value(0), "db_search_max_missed_cleavages")
			("db_search_max_dynamic_mods", boost::program_options::value<int>(&db_search_max_dynamic_mods)->default_value(0), "db_search_max_dynamic_mods")
			("db_search_min_enzymatic_termini", boost::program_options::value<int>(&db_search_min_enzymatic_termini)->default_value(2), "db_search_min_enzymatic_termini")
			("db_search_mass_tolerance", boost::program_options::value<double>(&db_search_mass_tolerance)->default_value(.05), "db_search_mass_tolerance")
			("db_search_PTM", boost::program_options::value<std::vector<DBPTM> >(&PTMs)->multitoken(), "PTMs")
			("db_search_enzyme", boost::program_options::value<std::vector<DBEnzyme>>(&enzymes)->multitoken(), "enzymes")
			("db_search_null_lambda", boost::program_options::value<double>(&null_lambda)->default_value(6), "db_search_null_lambda")
			("max_ms1_injection_time", boost::program_options::value<double>(&max_ms1_injection_time)->default_value(0.2), "max_ms1_injection_time")
			("max_ms2_injection_time", boost::program_options::value<double>(&max_ms2_injection_time)->default_value(0.5), "max_ms2_injection_time")
			("ms1_target_total_ion_count", boost::program_options::value<double>(&ms1_target_total_ion_count)->default_value(1e6), "ms1_target_total_ion_count")
			("ms2_target_total_ion_count", boost::program_options::value<double>(&ms2_target_total_ion_count)->default_value(1e5), "ms2_target_total_ion_count")
			;
	boost::program_options::variables_map vm_config;
	std::ifstream config_file(param_file_path.c_str());
	boost::program_options::store(boost::program_options::parse_config_file(config_file, config, true), vm_config);
	boost::program_options::notify(vm_config);
	//endregion


	ElutionShapeSimulator elution_shape_simulator(elution_tau, elution_sigma);
	std::unique_ptr<GroundTruthText> db = std::unique_ptr<GroundTruthText>(new GroundTruthText(sqlite_in_path, false));
	std::unique_ptr<Instrument> instrument = std::unique_ptr<Instrument>(new Instrument(resolution, dynamic_range, ms1_scan_time, ms2_scan_time,
																						scan_overhead_time, max_ms1_injection_time,
																						max_ms2_injection_time, ms1_target_total_ion_count,
																						ms2_target_total_ion_count));

	Sequencer sequencer(fasta_in_path, PTMs, enzymes, db_search_mass_tolerance, db_search_max_missed_cleavages, db_search_min_enzymatic_termini, db_search_min_mass, db_search_max_mass, db_search_max_dynamic_mods, null_lambda);
	Oracle oracle(db.get(), instrument.get(), elution_shape_simulator);
	MzMLWriter mzml_writer(mzml_out_path);
	FidoWriter fido_writer(fido_out_path);

	AcquisitionController* controller = get_controller(acquisition_algorithm_name, acquisition_param_values);

	auto start = std::chrono::high_resolution_clock::now();

	int ms1_count = 0;
	int ms2_count = 0;
	int quality_ms2_count = 0;
	double current_time = 0;

	std::cout << "Simulating Acquisition:" << std::endl;

	while (current_time <= acquisition_length) {
		std::unique_ptr<ScanRequest> scan_request = controller->get_scan_request(current_time);
		std::unique_ptr<Scan> scan = oracle.get_scan_data(scan_request.get(), current_time);

		if (scan->scan_type == Scan::ScanType::MS2) {
			ms2_count++;
			MS2Scan* tmp_scan = static_cast<MS2Scan*>(scan.get());
			sequencer.sequence_ms2_scan(tmp_scan);

			if (tmp_scan->probability >= 0 && tmp_scan->peptide != "DECOY") {
				fido_writer.write_peptide(tmp_scan->probability, tmp_scan->peptide, tmp_scan->proteins);
				if (tmp_scan->probability >= .9) quality_ms2_count++;
			}
		} else {
			ms1_count++;
		}

		controller->process_scan(scan.get());
		current_time += scan->elapsed_time;

		if (scan->scan_id % 20 == 0) {
			std::cout << "\rCurrent time: " << current_time << " seconds. MS1 count: " << ms1_count << ". MS2 count: " << ms2_count << ". Num PSMs >= 0.9: " << quality_ms2_count << std::flush;
		}

		mzml_writer.add_to_scan_buffer(std::move(scan));
		if (mzml_writer.buffer.size() > 100) mzml_writer.write_buffer();
	}

	std::cout << "\rCurrent time: " << current_time << " seconds. MS1 count: " << ms1_count << ". MS2 count: " << ms2_count << ". Num PSMs >= 0.9: " << quality_ms2_count << std::endl;

	mzml_writer.write_buffer();
	mzml_writer.output_file_end();

	mzml_writer.close_file();
	fido_writer.close_file();

	auto end = std::chrono::high_resolution_clock::now();

	std::cout << "Simulation Complete." << std::endl;
	std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << " seconds" << std::endl;

	return 0;
}