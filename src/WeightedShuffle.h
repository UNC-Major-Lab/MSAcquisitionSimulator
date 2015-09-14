//
// Created by Dennis Goldfarb on 8/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_WEIGHTEDSHUFFLE_H
#define MSACQUISITIONSIMULATOR_WEIGHTEDSHUFFLE_H

#include <algorithm>
#include <vector>
#include "Globals.h"

class WeightedShuffle {
private:
	static void update_normalized_prefix_sum(std::vector<double> &prefix_sum, int index);
	static void normalized_prefix_sum(std::vector<double> &weights);
	static void prefix_sum(std::vector<double> &weights);
	static void normalize_prefix_sum(std::vector<double> &prefix_sum);

public:
	template <class T> static std::vector<T> shuffle(std::vector<T> &values, std::vector<double> weights);
};




template<class T>
std::vector<T> WeightedShuffle::shuffle(std::vector<T> &values, std::vector<double> weights) {
	normalized_prefix_sum(weights);

	std::vector<T> shuffled;
	for (int i = 0; i < values.size(); ++i) {
		double random = g_uniform_distribution(g_rng);
		auto itr = std::upper_bound(weights.begin(), weights.end(), random);
		int random_index = (int) (itr-weights.begin());
		shuffled.push_back(values[random_index]);
		update_normalized_prefix_sum(weights, random_index);
	}
	return shuffled;
}


#endif //MSACQUISITIONSIMULATOR_WEIGHTEDSHUFFLE_H
