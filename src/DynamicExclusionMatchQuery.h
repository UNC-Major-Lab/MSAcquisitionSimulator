//
// Created by Dennis Goldfarb on 8/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_DYNAMICEXCLUSIONMATCHQUERY_H
#define MSACQUISITIONSIMULATOR_DYNAMICEXCLUSIONMATCHQUERY_H


class DynamicExclusionMatchQuery {
private:

public:

	DynamicExclusionMatchQuery(double mz, double mz_tolerance) : mz(mz), mz_tolerance(mz_tolerance) {}

	double mz;
	double mz_tolerance;
};


#endif //MSACQUISITIONSIMULATOR_DYNAMICEXCLUSIONMATCHQUERY_H
