//
// Created by Dennis Goldfarb on 6/29/15.
//

#include "Globals.h"

// TODO uncomment time seed
//boost::random::mt19937 g_rng((unsigned int) std::chrono::system_clock::now().time_since_epoch().count());

boost::random::mt19937 g_rng(0);
boost::random::uniform_real_distribution<double> g_uniform_distribution;

double MAX_MASS;
double MAX_MZ;
double MIN_MZ;
double PRUNE_THRESHOLD;