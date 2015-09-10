//
// Created by Dennis Goldfarb on 8/4/15.
//

#ifndef MSACQUISITIONSIMULATOR_IONIZATIONEFFICIENCYSIMULATOR_H
#define MSACQUISITIONSIMULATOR_IONIZATIONEFFICIENCYSIMULATOR_H


#include "PeptideDigested.h"
#include <boost/math/distributions/binomial.hpp>

class IonizationEfficiencySimulator {
public:

	IonizationEfficiencySimulator() {}

	double calc_ionization_efficiency(const Peptide &peptide);

	std::map<int,double> calc_charge_distribution(const Peptide &peptide);
};


#endif //MSACQUISITIONSIMULATOR_IONIZATIONEFFICIENCYSIMULATOR_H
