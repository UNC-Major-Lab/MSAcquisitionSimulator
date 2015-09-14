//
// Created by Dennis Goldfarb on 9/10/15.
//

#ifndef MSACQUISITIONSIMULATOR_IONDAO_H
#define MSACQUISITIONSIMULATOR_IONDAO_H

#include <string>

class IonDAO {
private:

public:
	IonDAO(int neutrons, int charge, int start, int end, double abundance, double rt,
		   double rt_start, double rt_end, double mz, std::string modified_sequence) :
			neutrons(neutrons), charge(charge), start(start), end(end), abundance(abundance),
			rt(rt), rt_start(rt_start), rt_end(rt_end), modified_sequence(modified_sequence),
			mz(mz), mass(mz*charge) {};

	int neutrons;
	int charge;
	int start;
	int end;
	double abundance;
	double rt;
	double rt_start;
	double rt_end;
	double mz;
	double mass;
	std::string modified_sequence;

	static bool less_rt_start(const IonDAO& a, const IonDAO& b);
};


#endif //MSACQUISITIONSIMULATOR_IONDAO_H
