//
// Created by Dennis Goldfarb on 8/18/15.
//

#ifndef MSACQUISITIONSIMULATOR_HISTOGRAM_H
#define MSACQUISITIONSIMULATOR_HISTOGRAM_H


#include <map>
#include <cmath>
#include <limits>
#include <iostream>

class Histogram {
private:

public:
	Histogram(std::string title, std::string y_axis, std::string x_axis) : title(title), y_axis(y_axis), x_axis(x_axis) {
		max_bin = std::numeric_limits<int>::min();
		min_bin = std::numeric_limits<int>::max();
	};
	~Histogram() {};

	std::string title;
	std::string y_axis;
	std::string x_axis;
	std::map<int,int> bin2count;
	int min_bin;
	int max_bin;

	void add_data(double d);
	void print_histogram();
};


#endif //MSACQUISITIONSIMULATOR_HISTOGRAM_H
