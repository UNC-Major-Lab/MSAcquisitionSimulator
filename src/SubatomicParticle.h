//
// Created by Dennis Goldfarb on 6/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_SUBATOMICPARTICLE_H
#define MSACQUISITIONSIMULATOR_SUBATOMICPARTICLE_H

#include <map>
#include <string>

struct SubatomicParticle {
	const double mass;
	const int charge;

	constexpr SubatomicParticle(double mass, int charge) : mass(mass), charge(charge) { }
};

namespace subatomic_particles {
	constexpr SubatomicParticle NEUTRON  {1.008664916,     0};
	constexpr SubatomicParticle PROTON   {1.007276466812,  1};
	constexpr SubatomicParticle ELECTRON {5.489e-4,       -1};

	const std::map<const std::string, const SubatomicParticle*> name2particle{{"n", &NEUTRON}, {"p", &PROTON},  {"e", &ELECTRON}};

}

#endif //MSACQUISITIONSIMULATOR_SUBATOMICPARTICLE_H
