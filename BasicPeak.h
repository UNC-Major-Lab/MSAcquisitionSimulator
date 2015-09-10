//
// Created by Dennis Goldfarb on 9/3/15.
//

#ifndef MSACQUISITIONSIMULATOR_BASICPEAK_H
#define MSACQUISITIONSIMULATOR_BASICPEAK_H


class BasicPeak {

private:

public:
	BasicPeak(double mz, double intensity) :
			mz(mz), intensity(intensity) {}

	double mz;
	double intensity;

	static bool greater_intensity(const BasicPeak& a, const BasicPeak& b);
	static bool less_mz(const BasicPeak& a, const BasicPeak& b);
};


#endif //MSACQUISITIONSIMULATOR_BASICPEAK_H
