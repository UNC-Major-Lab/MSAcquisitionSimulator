//
// Created by Dennis Goldfarb on 7/12/15.
//

#ifndef MSACQUISITIONSIMULATOR_RANDOMFACTORY_H
#define MSACQUISITIONSIMULATOR_RANDOMFACTORY_H


#include <chrono>
#include <boost/random.hpp>

class RandomFactory {

private:


public:

	RandomFactory();

	boost::random::mt19937 rng;
	boost::random::normal_distribution<double> normal_distribution;
	boost::random::uniform_real_distribution<double> uniform_distribution;

	double next_double_normal();
	double next_double_uniform();

};

#endif //MSACQUISITIONSIMULATOR_RANDOMFACTORY_H
