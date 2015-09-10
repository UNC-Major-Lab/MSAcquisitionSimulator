//
// Created by Dennis Goldfarb on 8/29/15.
//

#include "WeightedShuffle.h"


void WeightedShuffle::update_normalized_prefix_sum(std::vector<double> &prefix_sum, int index) {
	double loss = prefix_sum[index];
	if (index > 0) loss-= prefix_sum[index-1];

	for (; index < prefix_sum.size(); ++index) {
		prefix_sum[index] -= loss;
	}

	normalize_prefix_sum(prefix_sum);
}

void WeightedShuffle::normalized_prefix_sum(std::vector<double> &weights) {
	prefix_sum(weights);
	normalize_prefix_sum(weights);
}

void WeightedShuffle::prefix_sum(std::vector<double> &weights) {
	double prev = 0;
	for (int i = 0; i < weights.size(); i++) {
		prev += weights[i];
		weights[i] = prev;
	}
}

void WeightedShuffle::normalize_prefix_sum(std::vector<double> &prefix_sum) {
	if (prefix_sum[prefix_sum.size()-1] > 0) {
		double c_max = prefix_sum[prefix_sum.size()-1];
		for (int i = 0; i < prefix_sum.size(); i++) prefix_sum[i] /= c_max;
	}
}
