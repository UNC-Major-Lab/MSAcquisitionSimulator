#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <boost/program_options.hpp>
#include <boost/random.hpp>
#include <chrono>
#include "FastaParser.h"
#include "Protein.h"

#include "MSAcquisitionSimulatorConfig.h"

boost::random::mt19937 rng((unsigned int) std::chrono::system_clock::now().time_since_epoch().count());
boost::random::lognormal_distribution<double> lognormal_distribution;
boost::random::normal_distribution<double>  normal_distribution;


double get_random_abundance(double random_double) {
    double min = 1;

    if (random_double < min) {
        std::cerr << "Warning: random abundance was less than " << min << ". It is replaced with an abundance of " << min << std::endl;
        random_double = min;
    }
    return random_double;
}

std::vector<double> * sample_proteins(std::vector<Protein> &proteins, int sample_size, std::string distribution) {
	std::random_shuffle(proteins.begin(), proteins.end());
	std::vector<double>* abundances = new std::vector<double>((unsigned long) sample_size);

	if (distribution == "normal") {
		for (int i = 0; i < sample_size; ++i) (*abundances)[i] = get_random_abundance(normal_distribution(rng));
	} else if (distribution == "lognormal") {
		for (int i = 0; i < sample_size; ++i) (*abundances)[i] = get_random_abundance(lognormal_distribution(rng));
	}

	return abundances;
}

void print_stats(std::vector<double> &abundances) {
    std::sort(abundances.begin(), abundances.end());
    double min = *std::min_element(abundances.begin(), abundances.end());
    double max = *std::max_element(abundances.begin(), abundances.end());
    double dynamic_range = std::log10(max) - std::log10(min);

    double sum = std::accumulate(abundances.begin(), abundances.end(), 0.0);
    double mean = sum / abundances.size();
    double sq_sum = std::inner_product(abundances.begin(), abundances.end(), abundances.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / abundances.size() - mean * mean);

    double median = abundances[abundances.size()/2];

    std::cout << "Number of proteins sampled: " << abundances.size() << std::endl << std::endl;
    std::cout << "Abundance distribution statistics" << std::endl;
    std::cout << "Min: " << min << "\tMax: " << max << "\tDynamic range: " << dynamic_range << std::endl;
    std::cout << "Median: " << median << "\tMean: " << mean << "\tStdev: " << stdev << std::endl;
}

int main(int argc, const char ** argv) {

    std::cout << "MSAcquisitionSimulator version " << MSAcquisitionSimulator_VERSION_MAJOR << "." << MSAcquisitionSimulator_VERSION_MINOR << std::endl;
    std::cout << "FASTASampler version " << FASTASampler_VERSION_MAJOR << "." << FASTASampler_VERSION_MINOR << std::endl;


    double mean;
    double stdev;
    double percentage;
    int num_prot;
    std::string distribution;
    std::string fasta_out_path;
    std::string fasta_in_path;

    //region Command line specification
    boost::program_options::options_description general("USAGE: FASTASampler [options] input.fasta\n\nOptions");
    general.add_options()
            ("help", "Print usage and exit.")
            ("distribution,d", boost::program_options::value<std::string>(&distribution), "Choose abundance distribution. Options: normal, lognormal")
            ("mean,m", boost::program_options::value<double>(&mean)->default_value(10), "The mean of the normal distribution, or m parameter for log-normal (log10).")
            ("stdev,s", boost::program_options::value<double>(&stdev)->default_value(0.9), "The standard deviation of the normal distribution, or s parameter for log-normal (log10). 99.97% of the data will be within +-3 standard deviations")
            ("percentage,p", boost::program_options::value<double>(&percentage)->default_value(1), "Percentage of proteins to sample from. Takes precedence over --numprot.")
            ("numprot,n", boost::program_options::value<int>(&num_prot), "Number of proteins to sample from FASTA.")
            ("fasta_out,o", boost::program_options::value<std::string>(&fasta_out_path)->default_value("sample.sim.fasta"), "output path for sampled FASTA file.")
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

    if (distribution == "normal") {
        normal_distribution = boost::random::normal_distribution<double>(mean, stdev);
    } else if (distribution == "lognormal") {
        mean = std::log(std::pow(10,mean));
        stdev = std::log(std::pow(10,stdev));
        lognormal_distribution = boost::random::lognormal_distribution<double>(mean, stdev);
    } else {
        std::cout << "Invalid choice for distribution parameter" << std::endl;
        std::cout << general << std::endl;
        return 0;
    }
    //endregion


    std::vector<Protein> proteins = parse_FASTA(fasta_in_path.c_str());

    int sample_size = (int) proteins.size();
    std::cout << "Number of proteins in FASTA: " << sample_size << std::endl;
    if (vm.count("percentage")) {
        if (percentage <= 0 || percentage > 1) {
            std::cout << "Invalid number percentage: " << percentage << " must be  > 0 and <= 1" << std::endl;
            return 0;
        }
        sample_size = (int) round(proteins.size() * percentage);
    } else if (vm.count("numprot")) {
        sample_size = num_prot;
        if (sample_size < 0 || sample_size > proteins.size()) {
            std::cout << "Invalid number of proteins to sample from: " << sample_size << ". FASTA contains " << proteins.size() << " proteins." << std::endl;
            return 0;
        }
    }


	std::vector<double> *abundances = sample_proteins(proteins, sample_size, distribution);

	std::ofstream fasta_out(fasta_out_path.c_str());
	for (int i = 0; i < abundances->size(); i++) {
		fasta_out << ">" << proteins[i].name << " " << proteins[i].description << " #" << (*abundances)[i] << std::endl;
		fasta_out << proteins[i].sequence << std::endl;
	}

	fasta_out.close();
	print_stats(*abundances);

    return 0;
}