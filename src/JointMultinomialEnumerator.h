//
// Created by Dennis Goldfarb on 7/21/15.
//

#ifndef MSACQUISITIONSIMULATOR_JOINTMULTINOMIALENUMERATOR_H
#define MSACQUISITIONSIMULATOR_JOINTMULTINOMIALENUMERATOR_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include "PTM.h"

using namespace std;

struct KeyHasher {
	std::size_t operator()(const vector<short>& k) const
	{
		size_t my_hash = 0;
		for (int i = 0; i < k.size(); ++i) {
			my_hash ^= hash<short>()(k[i]) << (i%16);
		}
		return my_hash;
	}
};

template <class T>
class JointMultinomialEnumerator {

private:

	static bool comp2(const pair<vector<short>,double>& a, const pair<vector<short>,double>& b);
	static bool comp3(const vector<pair<double,T>>& a, const vector<pair<double,T>>& b);

public:

	JointMultinomialEnumerator(vector<vector<pair<double,T>>> probabilities, double threshold);
	~JointMultinomialEnumerator() {}

	int length_minus_one;
	double threshold;
	vector<vector<pair<double,T>>> log_probabilities;
	vector<double> cumulative_probabilities;
	vector<short> next_state;

	bool has_next();
	pair<vector<short>, double> next_combination();

};


template<class T>
JointMultinomialEnumerator<T>::JointMultinomialEnumerator(vector<vector<pair<double,T>>> probabilities, double threshold)
		: log_probabilities(probabilities), threshold(threshold) {
	sort(log_probabilities.begin(), log_probabilities.end(), comp3);

	length_minus_one = (int) log_probabilities.size()-1;

	if (log_probabilities.size() == 0) {
		cout << "List of states is empty" << endl;
		return;
	}

	cumulative_probabilities.push_back(0);
	for (int i=0; i<= length_minus_one; i++) {
		cumulative_probabilities.push_back(cumulative_probabilities[i]+log_probabilities[i][0].first);
	}

	if (std::exp(cumulative_probabilities[length_minus_one]) >= threshold) {
		next_state = vector<short>(log_probabilities.size(),0);
	} else {
		cout << "Nothing passes threshold" << endl;
	}
}

template<class T>
bool JointMultinomialEnumerator<T>::has_next() {
	return next_state.size() > 0;
}

template<class T>
bool JointMultinomialEnumerator<T>::comp3(const vector<pair<double,T>> &a, const vector<pair<double,T>> &b) {
	return a[0].first > b[0].first;
}

template<class T>
bool JointMultinomialEnumerator<T>::comp2(const pair<vector<short>, double> &a, const pair<vector<short>, double> &b) {
	for (int i = 0; i < a.first.size(); ++i) {
		if (a.first[i] > b.first[i]) return true;
		if (a.first[i] < b.first[i]) return false;
	}
	return true;
}

template <class T>
pair<vector<short>, double> JointMultinomialEnumerator<T>::next_combination() {
	pair<vector<short>, double> comb(next_state, std::exp(cumulative_probabilities[length_minus_one]));

	int i = length_minus_one;
	while (i <= length_minus_one && i >= 0) {
		if (next_state[i] < log_probabilities[i].size() - 1 || next_state[i] == -1) {
			next_state[i]++;
			cumulative_probabilities[i+1] = cumulative_probabilities[i] + log_probabilities[i][next_state[i]].first;

			if (std::exp(cumulative_probabilities[i+1]) < threshold) {
				next_state[i] = -1;
				i--;
			} else {
				i++;
			}
		} else {
			next_state[i] = -1;
			i--;
		}
	}
	if (i == -1) next_state.clear();

	return comb;
}

#endif //MSACQUISITIONSIMULATOR_JOINTMULTINOMIALENUMERATOR_H
