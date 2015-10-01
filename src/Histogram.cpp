//
// Created by Dennis Goldfarb on 8/18/15.
//

#include "Histogram.h"

void Histogram::add_data(double d) {
	if (d > 0 && !std::isinf(d)) {
		int bin = (int) floor(log10(d));
		if (bin2count.find(bin) == bin2count.end()) bin2count[bin] = 0;
		bin2count[bin]++;
	}
 }

void Histogram::print_histogram() {
	if (bin2count.size() == 0) return;
	std::cout << "\n" << title << std::endl;
	std::cout << "\t\t" << "x-axis: " << x_axis << std::endl;

	std::vector<int> keys;

	for (auto pair = bin2count.begin(); pair != bin2count.end(); ++pair) {
		keys.push_back(pair->first);
	}

	sort(keys.begin(), keys.end());

	for (auto itr = keys.begin(); itr != keys.end(); ++itr) {
		int bin = *itr;
		std::cout << bin << "\t";
		if (bin2count.find(bin) != bin2count.end()) {
			int max_count = (int) floor(log2(bin2count[bin]));
			for (int i = 0; i <= max_count; i++) {
				std::cout << "*";
			}
		}
		std::cout << std::endl;
	}
	std::cout << "y-axis: " << y_axis << std::endl;
}
