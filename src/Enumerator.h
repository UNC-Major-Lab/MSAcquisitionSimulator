//
// Created by Dennis Goldfarb on 9/6/15.
//

#ifndef MSACQUISITIONSIMULATOR_ENUMERATOR_H
#define MSACQUISITIONSIMULATOR_ENUMERATOR_H

#include <vector>

template <class T>
class Enumerator {
private:

public:
	Enumerator(std::vector<std::vector<T>> states, int max_states) : states(states), max_states(max_states){
		for (int i = 0; i < states.size(); i++) {
			current_state.push_back(0);
		}
		is_first = true;
	};

	std::vector<std::vector<T>> states;
	std::vector<short> current_state;

	int max_states;
	bool is_first;
	bool has_next();
	std::vector<short> next_combination();
};

template <class T>
bool Enumerator<T>::has_next() {
	for (int i = 0; i < current_state.size(); i++) {
		if (current_state[i] < states[i].size()-1) return true;
	}
	return is_first;
}

template <class T>
std::vector<short> Enumerator<T>::next_combination() {
	if (is_first) {
		is_first = false;
	} else if (current_state.size() > 0) {
		current_state[current_state.size() - 1]++;
		for (int i = (int) current_state.size() - 1; i >= 0; i--) {
			if (current_state[i] == states[i].size()) {
				current_state[i] = 0;
				current_state[i - 1]++;
			} else {
				return current_state;
			}
		}
	}
	return current_state;
}


#endif //MSACQUISITIONSIMULATOR_ENUMERATOR_H
