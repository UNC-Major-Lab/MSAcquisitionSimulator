//
// Created by Dennis Goldfarb on 6/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_GLOBALS_H
#define MSACQUISITIONSIMULATOR_GLOBALS_H


#include <iostream>
#include <map>
#include <boost/random.hpp>

enum class TERMINUS {N_TERM, C_TERM};

enum class NEUTRAL_LOSSES {NH3, H2O};

// TODO uncomment time seed
//extern boost::random::mt19937 g_rng((unsigned int) std::chrono::system_clock::now().time_since_epoch().count());
extern boost::random::mt19937 g_rng;
extern boost::random::uniform_real_distribution<double> g_uniform_distribution;

extern double PRUNE_THRESHOLD;
extern double MIN_MZ;
extern double MAX_MZ;
extern double MAX_MASS;



#endif //MSACQUISITIONSIMULATOR_GLOBALS_H
