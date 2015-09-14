//
// Created by Dennis Goldfarb on 6/28/15.
//
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <chrono>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/bind.hpp>

#include "BioLCCC/biolccc.h"

#include "DefaultPTM.h"
#include "DefaultEnzyme.h"
#include "FASTAParser.h"
#include "PeptideGenerator.h"
#include "IonizationSimulator.h"
#include "Histogram.h"
#include "GroundTruthText.h"

#include "MSAcquisitionSimulatorConfig.h"

/*
 * Begin: Command line and config file parsing for user defined classes
 */
void validate(boost::any& v, const std::vector<std::string>& values, TERMINUS* target_type, int) {
	if (values[0][0] == 'N') v = TERMINUS::N_TERM;
	else if (values[0][0] == 'C') v = TERMINUS::C_TERM;
}

void validate(boost::any& v, const std::vector<std::string>& values, DefaultEnzyme* target_type, int) {
	double cleavage_probability;
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
			("cleavage_probability,p", boost::program_options::value<double>(&cleavage_probability), "")
			("blocking_residues,b", boost::program_options::value<std::string>(&blocking_residues), "")
			("terminus,t", boost::program_options::value<TERMINUS>(&terminus), "")
			;

	boost::tokenizer<separator_type> tokens(values[0], separator);
	std::vector<std::string> result;

	std::copy(tokens.begin(), tokens.end(), std::back_inserter(result));

	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::command_line_parser(result).options(general).run(), vm);
	boost::program_options::notify(vm);
	std::cout << name << " " << residues << " " << cleavage_probability << " " << blocking_residues  << std::endl;
	v = DefaultEnzyme(name, residues, blocking_residues, terminus, cleavage_probability);
}

void validate(boost::any& v, const std::vector<std::string>& values, DefaultPTM* target_type, int) {
	double site_probability;
	double abundance;
	double bind_energy;
	bool blocks_cleavage;
	bool post_digestion;
	std::string name;
	std::string abbreviation;
	std::string residues;
	std::string formula;

	typedef boost::escaped_list_separator<char> separator_type;
	separator_type separator("\\",    // The escape characters.
							 "\t ",    // The separator characters.
							 "\"\'"); // The quote characters.

	boost::program_options::options_description general("PTMs");
	general.add_options()
			("name,n", boost::program_options::value<std::string>(&name), "")
			("abbreviation,b", boost::program_options::value<std::string>(&abbreviation), "")
			("residues,r", boost::program_options::value<std::string>(&residues), "")
			("site_probabiltiy,p", boost::program_options::value<double>(&site_probability), "")
			("abundance,a", boost::program_options::value<double>(&abundance), "")
			("formula,m", boost::program_options::value<std::string>(&formula), "" )
			("bind_energy,e", boost::program_options::value<double>(&bind_energy), "")
			("blocks_cleavage", boost::program_options::bool_switch(&blocks_cleavage),  "")
			("post_digestion", boost::program_options::bool_switch(&post_digestion),"")
			;

	boost::tokenizer<separator_type> tokens(values[0], separator);
	std::vector<std::string> result;

	std::copy(tokens.begin(), tokens.end(), std::back_inserter(result));

	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::command_line_parser(result).options(general).run(), vm);
	boost::program_options::notify(vm);

	double bind_area = (bind_energy == 0) ? 0 : 1;

	v = DefaultPTM(name, abbreviation, residues, formula, site_probability, bind_energy, bind_area, abundance, blocks_cleavage, post_digestion);
}
/*
 * End: Command line and config file parsing for user defined classes
 */



/*
 * PTMs need to be added to BioLCCC for retention time prediction
 */
void register_ptm(PTM *ptm, BioLCCC::ChemicalBasis &chemical_basis) {
	for (char residue : ptm->residues) {
		std::string pep_abbr = ptm->abbreviation;
		if (residue != 'n' && residue != 'c') pep_abbr += residue;
		chemical_basis.addChemicalGroup(BioLCCC::ChemicalGroup(ptm->name, pep_abbr, ptm->bind_energy,
															   ptm->molecular_formula.get_monoisotopic_mass(),
															   ptm->molecular_formula.get_average_mass(),
															   ptm->bind_area));
	}
}


/*
 * Main function to perform all steps of the simulation
 */
void generate_ground_truth(IonizationEfficiencySimulator &ionization_efficiency_simulator, ElutionShapeSimulator &elution_shape_simulator,
						   std::string fasta_in_path, PeptideGenerator &peptide_generator, BioLCCC::ChromoConditions &chromo_conditions,
						   BioLCCC::ChemicalBasis &chemical_basis, double gradient_duration, GroundTruthText &db) {

	std::vector<Protein> proteins = parse_FASTA(fasta_in_path.c_str());
	std::map<Peptide,double> peptides;
	std::map<Protein*, double> protein2max_abundance;

	std::cout << std::endl;

	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < proteins.size(); ++i) {
		std::cout << "\rNumber of proteins processed: " << i << " of " << proteins.size() << ". Currently processing: "
		<< proteins[i].name << ". Abundance: " << proteins[i].abundance << ". Sequence length: " << proteins[i].sequence.length() << std::flush;

		protein2max_abundance[&proteins[i]] = 0;

		std::vector<Peptide> peps = peptide_generator.generate_peptides(proteins[i]);
		for (Peptide &p: peps) {
			if (peptides.find(p) == peptides.end()) peptides[p] = p.abundance;
			else peptides.at(p)+= p.abundance;
		}
	}

	std::cout << "\rNumber of proteins processed: " << proteins.size() << " of " << proteins.size() << std::endl;

	Histogram ion_histogram("Ion abundance distribution", "log10(ion abundance)", "log2(count)");

	int num_processed_peptides = 0;
	int num_processed_ions = 0;

	for (std::pair<const Peptide, double> &pair : peptides) {
		if (num_processed_peptides%1000 == 0) {
			std::cout << "\rNumber of peptides processed: " << num_processed_peptides << " of " << peptides.size() << ". Number of ions passing abundance thresholds: " << num_processed_ions << std::flush;
		}
		++num_processed_peptides;

		pair.second *= ionization_efficiency_simulator.calc_ionization_efficiency(pair.first);
		if (pair.second <= PRUNE_THRESHOLD) continue;

		double RT_fast = BioLCCC::calculateRT(pair.first.get_modified_sequence(), chemical_basis, chromo_conditions, 21);

		if (RT_fast > gradient_duration+2) continue; // allow the RT center to be at most 2 minutes after the end of the gradient

		std::map<int,double> charge_to_percentage = ionization_efficiency_simulator.calc_charge_distribution(pair.first);

		for (std::pair<int,double> pair2 : charge_to_percentage) {
			double new_abundance = pair2.second * pair.second;
			if (new_abundance <= PRUNE_THRESHOLD) continue;

			protein2max_abundance[pair.first.protein] = std::max(protein2max_abundance[pair.first.protein], new_abundance);

			ion_histogram.add_data(new_abundance);

			std::vector<double> isotope_mz;
			std::vector<double> isotope_abundance;
			mercury::mercury(isotope_mz, isotope_abundance, pair.first.get_composition(pair2.first), pair2.first, 10e-30);

			db.insert_ions(new_abundance, pair2.first, RT_fast*60, &pair.first, isotope_mz, isotope_abundance, elution_shape_simulator);

			++num_processed_ions;
		}
	}
	std::cout << "\rNumber of peptides processed: " << peptides.size() << " of " << peptides.size() << ". Number of ions passing abundance thresholds: " << num_processed_ions << std::endl;

	auto end = std::chrono::high_resolution_clock::now();

	std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << " seconds" << std::endl;

	Histogram max_int_histogram("Most abundant ion for protein", "log10(ion abundance)", "log2(count)");
	for (std::pair<Protein*,double> pair : protein2max_abundance) {
		max_int_histogram.add_data(pair.second);
	}

	ion_histogram.print_histogram();
	max_int_histogram.print_histogram();
}



int main(int argc, const char ** argv) {

	std::cout << "MSAcquisitionSimulator version " << MSAcquisitionSimulator_VERSION_MAJOR << "." << MSAcquisitionSimulator_VERSION_MINOR << "." << MSAcquisitionSimulator_VERSION_PATCH << std::endl;
	std::cout << "GroundTruthSimulator version " << GroundTruthSimulator_VERSION_MAJOR << "." << GroundTruthSimulator_VERSION_MINOR << "." << GroundTruthSimulator_VERSION_PATCH << std::endl;


	std::string param_file_path;
	std::string sqlite_out_path;
	std::string fasta_in_path;

	double column_length;
	double column_diameter;
	double column_pore_size;
	double second_solvent_concentration_a;
	double second_solvent_concentration_b;
	double gradient_percent_b_start;
	double gradient_percent_b_end;
	double gradient_duration;
	double gradient_flow_rate;

	double elution_tau;
	double elution_sigma;

	std::vector<DefaultPTM> PTMs;
	std::vector<DefaultEnzyme> enzymes;

	BioLCCC::ChromoConditions chromo_conditions;
	BioLCCC::ChemicalBasis chemical_basis(BioLCCC::RP_ACN_FA_ROD);

	//region Command line specification
	boost::program_options::options_description general("USAGE: GroundTruthSimulator [options] input.fasta\n\nOptions");
	general.add_options()
			("help", "Print usage and exit.")
			("config,c", boost::program_options::value<std::string>(&param_file_path)->default_value("ground_truth.conf"), "Input path to config file.")
			("sqlite_out_path,o", boost::program_options::value<std::string>(&sqlite_out_path)->default_value("sample.sqlite"), "output path for sampled FASTA file.")
			;

	boost::program_options::options_description hidden("");
	hidden.add_options()
			("fasta_in", boost::program_options::value<std::string>(&fasta_in_path), "input path for FASTA file to sample from.")
			;

	boost::program_options::options_description all("Allowed options");
	all.add(general).add(hidden);

	boost::program_options::positional_options_description p;
	p.add("fasta_in", -1);

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
			("PTM", boost::program_options::value<std::vector<DefaultPTM> >(&PTMs)->multitoken(), "PTMs")
			("Enzyme", boost::program_options::value<std::vector<DefaultEnzyme> >(&enzymes)->multitoken(), "Enzymes")
			("column_length", boost::program_options::value<double>(&column_length), "column_length")
			("column_diameter", boost::program_options::value<double>(&column_diameter), "column_diameter")
			("column_pore_size", boost::program_options::value<double>(&column_pore_size), "column_pore_size")
			("second_solvent_concentration_a", boost::program_options::value<double>(&second_solvent_concentration_a), "second_solvent_concentration_a")
			("second_solvent_concentration_b", boost::program_options::value<double>(&second_solvent_concentration_b), "second_solvent_concentration_b")
			("gradient_percent_b_start", boost::program_options::value<double>(&gradient_percent_b_start), "gradient_percent_b_start")
			("gradient_percent_b_end", boost::program_options::value<double>(&gradient_percent_b_end), "gradient_percent_b_end")
			("gradient_duration", boost::program_options::value<double>(&gradient_duration), "gradient_duration")
			("gradient_flow_rate", boost::program_options::value<double>(&gradient_flow_rate), "gradient_flow_rate")
			("max_mass", boost::program_options::value<double>(&MAX_MASS)->default_value(10000), "max_mass")
			("max_mz", boost::program_options::value<double>(&MAX_MZ)->default_value(3000), "max_mz")
			("min_mz", boost::program_options::value<double>(&MIN_MZ)->default_value(200), "min_mz")
			("prune_threshold", boost::program_options::value<double>(&PRUNE_THRESHOLD)->default_value(1e4), "prune_threshold")
			("elution_tau", boost::program_options::value<double>(&elution_tau)->default_value(4), "elution_tau")
			("elution_sigma", boost::program_options::value<double>(&elution_sigma)->default_value(6), "elution_sigma")
			;
	boost::program_options::variables_map vm_config;
	std::ifstream config_file(param_file_path.c_str());
	boost::program_options::store(boost::program_options::parse_config_file(config_file, config, true), vm_config);
	boost::program_options::notify(vm_config);
	//endregion

	chromo_conditions.setColumnLength(column_length);
	chromo_conditions.setColumnDiameter(column_diameter);
	chromo_conditions.setColumnPoreSize(column_pore_size);
	chromo_conditions.setSecondSolventConcentrationA(second_solvent_concentration_a);
	chromo_conditions.setSecondSolventConcentrationB(second_solvent_concentration_b);
	chromo_conditions.setGradient(BioLCCC::Gradient(gradient_percent_b_start, gradient_percent_b_end, gradient_duration));
	chromo_conditions.setFlowRate(gradient_flow_rate);

	GroundTruthText db(sqlite_out_path, true);
	PeptideGenerator peptide_generator = PeptideGenerator(PTMs, enzymes);
	IonizationEfficiencySimulator ionization_efficiency_simulator;
	ElutionShapeSimulator elution_shape_simulator(elution_tau, elution_sigma);

	for (PTM* ptm : peptide_generator.modification_simulator.pre_digestion_modifications) {
		register_ptm(ptm, chemical_basis);
	}
	for (PTM* ptm : peptide_generator.modification_simulator.post_digestion_modifications) {
		register_ptm(ptm, chemical_basis);
	}

	generate_ground_truth(ionization_efficiency_simulator, elution_shape_simulator, fasta_in_path,
						  peptide_generator, chromo_conditions, chemical_basis, gradient_duration, db);


	std::cout << "Sorting ions by retention time. This might take a while.." << std::endl;
	db.write_sorted_file(); // Must be sorted by retention time to make the AcquisitionSimulator faster.
	return 0;
}