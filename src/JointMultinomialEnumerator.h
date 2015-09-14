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
	double calc_probability(vector<short> &state);
	void next_steps(const vector<short> &state, vector<vector<short>> &next_states);

public:

	JointMultinomialEnumerator(vector<vector<pair<double,T>>> probabilities, double threshold);
	~JointMultinomialEnumerator() {}

	int length_minus_one;
	bool multi = false;
	double threshold;
	vector<vector<pair<double,T>>> log_probabilities;
	unordered_set<vector<short>, KeyHasher> current_states_hash;
	vector<pair<vector<short>,double>> min_heap;
	vector<vector<short>> new_states;

	bool has_next();
	pair<vector<short>, double> next_combination();

};


template<class T>
JointMultinomialEnumerator<T>::JointMultinomialEnumerator(vector<vector<pair<double,T>>> probabilities, double threshold)
		: log_probabilities(probabilities), threshold(threshold) {
	for (vector<pair<double,T>> &v : log_probabilities) {
		if (v.size() > 2) multi = true;
	}
	sort(log_probabilities.begin(), log_probabilities.end(), comp3);

	/*for (const vector<double> &v : log_probabilities) {
		for (const double &p : v) {
			cout << exp(p) << " ";
		}
		cout << endl;
	}*/
	length_minus_one = (int) log_probabilities.size()-1;

	vector<short> init(log_probabilities.size(),0);
	if (log_probabilities.size() == 0) {
		cout << "List of states is empty" << endl;
		return;
	}
	double prob = calc_probability(init);
	if (prob >= threshold) {
		current_states_hash.insert(init);
		min_heap.push_back({init, prob});
	} else {
		cout << "Nothing passes threshold" << endl;
	}
}

template<class T>
bool JointMultinomialEnumerator<T>::has_next() {
	return min_heap.size() > 0;
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

template<class T>
double JointMultinomialEnumerator<T>::calc_probability(vector<short> &state) {
	double prob = 0;
	for (int i = 0; i < state.size(); ++i) {
		prob += log_probabilities[i][state[i]].first;
	}
	return exp(prob);
}

template <class T>
pair<vector<short>, double> JointMultinomialEnumerator<T>::next_combination() {
	vector<short> state = min_heap[0].first;

	pair<vector<short>, double> comb(min_heap[0]);

	current_states_hash.erase(state);
	pop_heap(min_heap.begin(), min_heap.end(), comp2); min_heap.pop_back();

	next_steps(state, new_states);

	for (vector<short> &ns : new_states) {
		if (current_states_hash.find(ns) == current_states_hash.end()) {
			double prob = calc_probability(ns);
			if (prob >= threshold) {
				current_states_hash.insert(ns);
				min_heap.push_back({ns, prob});
				push_heap(min_heap.begin(), min_heap.end(), comp2);
			}
		}
	}

	return comb;
}

template<class T>
void JointMultinomialEnumerator<T>::next_steps(const vector<short> &state, vector<vector<short>> &next_states) {
	next_states.clear();

	if (state[length_minus_one] == 0) {
		vector<short> next_state(state);
		next_state[length_minus_one] = 1;
		next_states.push_back(next_state);
	}

	for (int i = length_minus_one; i>0; --i) {
		if (state[i] == 1 && state[i-1] == 0) {
			vector<short> next_state(state);
			next_state[i] = 0;
			next_state[i-1] = 1;
			next_states.push_back(next_state);
			break;
		}
	}

	if (multi) {
		for (int i = length_minus_one; i >= 0; --i) {
			if (state[i] > 0 && state[i] < log_probabilities[i].size() - 1) {
				vector<short> next_state(state);
				next_state[i]++;
				next_states.push_back(next_state);
			}
		}
	}
}

#endif //MSACQUISITIONSIMULATOR_JOINTMULTINOMIALENUMERATOR_H
