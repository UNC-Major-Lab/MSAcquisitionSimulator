//
// Created by Dennis Goldfarb on 8/4/15.
//

#include "IonizationSimulator.h"
#include "Globals.h"

double IonizationEfficiencySimulator::calc_ionization_efficiency(const Peptide &peptide) {
	return g_uniform_distribution(g_rng);
}

std::map<int, double> IonizationEfficiencySimulator::calc_charge_distribution(const Peptide &peptide) {
	std::map<int,double> charge_to_percentage;


	int num_basic = 1+peptide.num_basic_residues();
	boost::math::binomial bd(num_basic, .7+(g_uniform_distribution(g_rng)*.3)); // success rate in [0.7,1.0]

	for (int i = 1; i <= num_basic; i++) {
		charge_to_percentage[i] = boost::math::pdf(bd,i);
		//std::cout << i << " " << num_basic << " " << boost::math::pdf(bd,i) << std::endl;
	}

	return charge_to_percentage;
}
