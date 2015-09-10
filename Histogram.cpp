//
// Created by Dennis Goldfarb on 8/18/15.
//

#include "Histogram.h"

void Histogram::add_data(double d) {
	int bin = (int) floor(log10(d));
	if (bin2count.find(bin) == bin2count.end()) bin2count[bin] = 0;
	max_bin = std::max(bin, max_bin);
	min_bin = std::min(bin, min_bin);
	bin2count[bin]++;
 }

void Histogram::print_histogram() {

	std::cout << "\n" << title << std::endl;
	std::cout << "\t\t" << "x-axis: " << x_axis << std::endl;
	for (int bin = max_bin; bin >= min_bin; bin--) {
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
