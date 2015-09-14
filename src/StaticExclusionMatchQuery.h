//
// Created by Dennis Goldfarb on 8/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_STATICEXCLUSIONMATCHQUERY_H
#define MSACQUISITIONSIMULATOR_STATICEXCLUSIONMATCHQUERY_H


class StaticExclusionMatchQuery {
private:

public:

	StaticExclusionMatchQuery(double mz, double mz_tolerance, double time) :
			mz(mz), mz_tolerance(mz_tolerance), time(time) {}

	double mz;
	double mz_tolerance;
	double time;
};


#endif //MSACQUISITIONSIMULATOR_STATICEXCLUSIONMATCHQUERY_H
