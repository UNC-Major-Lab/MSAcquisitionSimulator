//
// Created by Dennis Goldfarb on 7/12/15.
//

#include "RandomFactory.h"

RandomFactory::RandomFactory() {
	rng = boost::random::mt19937((unsigned int) std::chrono::system_clock::now().time_since_epoch().count());
	normal_distribution = boost::random::normal_distribution<double>(0,1);
};

double RandomFactory::next_double_uniform() {
	return uniform_distribution(rng);
}

double RandomFactory::next_double_normal() {
	return normal_distribution(rng);
}
