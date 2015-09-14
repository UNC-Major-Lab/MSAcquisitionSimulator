//
// Created by Dennis Goldfarb on 7/23/15.
//

#include <math.h>
#include "PeptideDigested.h"

double PeptideDigested::calc_digestion_probability() {
	return exp(log_C_cleavage_prob + log_N_cleavage_prob + log_no_cleavage_prob);
}
