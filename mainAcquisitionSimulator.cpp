//
// Created by Dennis Goldfarb on 8/19/15.
//

#include <iostream>
#include <fstream>
#include <chrono>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/bind.hpp>
#include "sqlite3.h"
#include "Globals.h"
#include "AcquisitionController.h"
#include "Oracle.h"
#include "GroundTruthDatabase.h"
#include "Sequencer.h"
#include "TopN.h"
#include "MzMLWriter.h"
#include "DBEnzyme.h"
#include "DBPTM.h"
#include "FidoWriter.h"

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

	std::cout << name << " " << residues << " " << blocking_residues  << std::endl;

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

	std::cout << name << " " << residues << " " << formula << std::endl;

	v = DBPTM(name, abbreviation, residues, formula, is_static);
}

int main(int argc, const char ** argv) {

	std::string sqlite_in_path;
	std::string mzml_out_path;
	std::string param_file_path;
	std::string fido_out_path;

	std::string acquisition_algorithm_name;
	std::string fasta_in_path;
	int num_ms2;
	int ms1_scan_time;
	int ms2_scan_time;
	int acquisition_length;
	double ms1_min_mz;
	double ms1_max_mz;
	double ms2_isolation_width;
	double dynamic_exclusion_tolerance;
	int dynamic_exclusion_time;
	bool dynamic_exclusion_enabled;

	double elution_tau;
	double elution_sigma;

	double max_ms1_injection_time;
	double max_ms2_injection_time;
	double ms1_target_total_ion_count;
	double ms2_target_total_ion_count;
	double resolution;
	double dynamic_range;

	double db_search_min_mass;
	double db_search_max_mass;
	double db_search_mass_tolerance;
	int db_search_max_missed_cleavages;
	int db_search_max_dynamic_mods;
	int db_search_min_enzymatic_termini;
	double null_lambda;

	std::vector<DBPTM> PTMs;
	std::vector<DBEnzyme> enzymes;


	//region Command line specification
	boost::program_options::options_description general("USAGE: AcquisitionSimulator [options] ground_truth.sqlite\n\nOptions");
	general.add_options()
			("help", "Print usage and exit.")
			("param,p", boost::program_options::value<std::string>(&param_file_path)->default_value("ground_truth.params"), "Input path to parameter file.")
			("mzml_out_path,o", boost::program_options::value<std::string>(&mzml_out_path)->default_value("sample.mzML"), "output path for mzML file.")
			("fido_out_path,f", boost::program_options::value<std::string>(&fido_out_path)->default_value("sample.fido"), "output path for fido file.")
			;

	boost::program_options::options_description hidden("");
	hidden.add_options()
			("sqlite_in_path", boost::program_options::value<std::string>(&sqlite_in_path), "input path for ground truth SQLite file made by GroundTruthSimulator.")
			;

	boost::program_options::options_description all("Allowed options");
	all.add(general).add(hidden);

	boost::program_options::positional_options_description p;
	p.add("sqlite_in_path", -1);

	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(all).positional(p).run(), vm);
	boost::program_options::notify(vm);
	//endregion

	//region config file specification
	boost::program_options::options_description config("Configuration file options");
	config.add_options()
			("acquisition_algorithm", boost::program_options::value<std::string>(&acquisition_algorithm_name)->default_value("TopN"), "acquisition algorithm")
			("num_ms2", boost::program_options::value<int>(&num_ms2)->default_value(10), "num_ms2")
			("ms1_scan_time", boost::program_options::value<int>(&ms1_scan_time)->default_value(128), "ms1_scan_time")
			("ms2_scan_time", boost::program_options::value<int>(&ms2_scan_time)->default_value(64), "ms2_scan_time")
			("ms1_min_mz", boost::program_options::value<double>(&ms1_min_mz)->default_value(200), "ms1_min_mz")
			("ms1_max_mz", boost::program_options::value<double>(&ms1_max_mz)->default_value(3000), "ms1_max_mz")
			("ms2_isolation_width", boost::program_options::value<double>(&ms2_isolation_width)->default_value(1), "ms2_isolation_width")
			("dynamic_exclusion_enabled", boost::program_options::value<bool>(&dynamic_exclusion_enabled)->default_value(true), "dynamic_exclusion_enabled")
			("dynamic_exclusion_tolerance", boost::program_options::value<double>(&dynamic_exclusion_tolerance)->default_value(0.05), "dynamic_exclusion_tolerance")
			("dynamic_exclusion_time", boost::program_options::value<int>(&dynamic_exclusion_time)->default_value(30), "dynamic_exclusion_time")
			("acquisition_length", boost::program_options::value<int>(&acquisition_length)->default_value(3600), "acquisition_length")
			("fasta", boost::program_options::value<std::string>(&fasta_in_path), "acquisition fasta_in_path")
			("elution_tau", boost::program_options::value<double>(&elution_tau)->default_value(4), "elution_tau")
			("elution_sigma", boost::program_options::value<double>(&elution_sigma)->default_value(6), "elution_sigma")
			("max_ms1_injection_time", boost::program_options::value<double>(&max_ms1_injection_time)->default_value(0.2), "max_ms1_injection_time")
			("max_ms2_injection_time", boost::program_options::value<double>(&max_ms2_injection_time)->default_value(0.5), "max_ms2_injection_time")
			("ms1_target_total_ion_count", boost::program_options::value<double>(&ms1_target_total_ion_count)->default_value(1e6), "ms1_target_total_ion_count")
			("ms2_target_total_ion_count", boost::program_options::value<double>(&ms2_target_total_ion_count)->default_value(1e5), "ms2_target_total_ion_count")
			("resolution", boost::program_options::value<double>(&resolution)->default_value(60000), "resolution")
			("dynamic_range", boost::program_options::value<double>(&dynamic_range)->default_value(5000), "dynamic_range")
			("db_search_min_mass", boost::program_options::value<double>(&db_search_min_mass)->default_value(200), "db_search_min_mass")
			("db_search_max_mass", boost::program_options::value<double>(&db_search_max_mass)->default_value(9000), "db_search_max_mass")
			("db_search_max_missed_cleavages", boost::program_options::value<int>(&db_search_max_missed_cleavages)->default_value(0), "db_search_max_missed_cleavages")
			("db_search_max_dynamic_mods", boost::program_options::value<int>(&db_search_max_dynamic_mods)->default_value(0), "db_search_max_dynamic_mods")
			("db_search_min_enzymatic_termini", boost::program_options::value<int>(&db_search_min_enzymatic_termini)->default_value(2), "db_search_min_enzymatic_termini")
			("db_search_mass_tolerance", boost::program_options::value<double>(&db_search_mass_tolerance)->default_value(.05), "db_search_mass_tolerance")
			("db_search_PTM", boost::program_options::value<std::vector<DBPTM> >(&PTMs)->multitoken(), "PTMs")
			("db_search_enzyme", boost::program_options::value<std::vector<DBEnzyme> >(&enzymes)->multitoken(), "enzymes")
			("null_lambda", boost::program_options::value<double>(&null_lambda)->default_value(6), "null_lambda")
			;
	boost::program_options::variables_map vm_config;
	std::ifstream config_file(param_file_path.c_str());
	boost::program_options::store(boost::program_options::parse_config_file(config_file, config, true), vm_config);
	boost::program_options::notify(vm_config);
	//endregion


	ElutionShapeSimulator elution_shape_simulator(elution_tau, elution_sigma);
	std::unique_ptr<GroundTruthDatabase> db = std::unique_ptr<GroundTruthDatabase>(new GroundTruthDatabase(sqlite_in_path, false));
	std::unique_ptr<Instrument> instrument = std::unique_ptr<Instrument>(new Instrument(resolution, dynamic_range));
	TopNParameters params(ms2_isolation_width, ms1_min_mz, ms1_max_mz, dynamic_exclusion_tolerance,
						  dynamic_exclusion_time, num_ms2, dynamic_exclusion_enabled, max_ms1_injection_time, max_ms2_injection_time,
						  ms1_target_total_ion_count, ms2_target_total_ion_count);
	std::unique_ptr<AcquisitionController> controller = std::unique_ptr<AcquisitionController>(new TopN(params));
	Sequencer sequencer(fasta_in_path, PTMs, enzymes, db_search_mass_tolerance, db_search_max_missed_cleavages, db_search_min_enzymatic_termini, db_search_min_mass, db_search_max_mass, db_search_max_dynamic_mods, null_lambda);
	Oracle oracle(db.get(), instrument.get(), elution_shape_simulator);
	MzMLWriter writer(mzml_out_path);
	FidoWriter fido_writer(fido_out_path);

	auto start = std::chrono::high_resolution_clock::now();

	double current_time = 3600;

	//while (current_time <= acquisition_length) {
	while (current_time <= 3700.0) {
		std::unique_ptr<ScanRequest> scan_request = controller->get_scan_request(current_time);
		std::unique_ptr<Scan> scan = oracle.get_scan_data(scan_request.get(), current_time);

		if (scan->scan_type == Scan::ScanType::MS2) {
			MS2Scan* tmp_scan = static_cast<MS2Scan*>(scan.get());
			sequencer.sequence_ms2_scan(tmp_scan);

			if (tmp_scan->probability > 0 && tmp_scan->peptide != "DECOY") {
				fido_writer.write_peptide(tmp_scan->probability, tmp_scan->peptide, tmp_scan->proteins);
			}
		} else {
			std::cout << "MS1!" << std::endl;
		}

		controller->process_scan(scan.get());
		current_time += scan->elapsed_time;

		if (scan->scan_id % 100 == 0) {
			std::cout << "Current time: " << current_time << " Elapsed time: " << scan->elapsed_time << " ScanType: " << scan->scan_type << " Num peaks: " << scan->peaks.size() << std::endl;
		}

		writer.add_to_write_buffer(std::move(scan));
		if (writer.buffer.size() > 100) writer.write_buffer();
	}

	writer.write_buffer();
	writer.output_file_end();
	writer.close_file();
	fido_writer.close_file();

	auto end = std::chrono::high_resolution_clock::now();

	std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << " seconds" << std::endl;


	return 0;
}