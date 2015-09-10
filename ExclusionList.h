//
// Created by Dennis Goldfarb on 8/28/15.
//

#ifndef MSACQUISITIONSIMULATOR_EXCLUSIONLIST_H
#define MSACQUISITIONSIMULATOR_EXCLUSIONLIST_H

#include <set>

template <class T>
class ExclusionList {
private:

public:
	ExclusionList(double mz_tolerance) : mz_tolerance(mz_tolerance) {}

	void add_entry(T entry);

	double mz_tolerance;
	std::multiset<T> entries;
};


template <class T>
void ExclusionList<T>::add_entry(T entry) {
	entries.insert(entry);
}

#endif //MSACQUISITIONSIMULATOR_EXCLUSIONLIST_H
