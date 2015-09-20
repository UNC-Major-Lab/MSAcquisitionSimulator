//
// Created by Dennis Goldfarb on 8/18/15.
//

#ifndef MSACQUISITIONSIMULATOR_HISTOGRAM_H
#define MSACQUISITIONSIMULATOR_HISTOGRAM_H


#include <map>
#include <cmath>
#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>

class Histogram {
private:

public:
	Histogram(std::string title, std::string y_axis, std::string x_axis) : title(title), y_axis(y_axis), x_axis(x_axis) {
	};
	~Histogram() {};

	std::string title;
	std::string y_axis;
	std::string x_axis;
	std::map<int,int> bin2count;

	void add_data(double d);
	void print_histogram();
};


#endif //MSACQUISITIONSIMULATOR_HISTOGRAM_H
