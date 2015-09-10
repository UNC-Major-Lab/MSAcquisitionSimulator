//
// Created by Dennis Goldfarb on 6/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_RESIDUE_H
#define MSACQUISITIONSIMULATOR_RESIDUE_H

#include <map>
#include <string>
#include "Molecule.h"

class Residue : public Molecule {

private:


public:
	double monoisotopic_mass;

	Residue() {};

	Residue(double monoisotopic_mass, std::string formula) : Molecule(formula), monoisotopic_mass(monoisotopic_mass) {};

};

namespace residues {

	extern Residue A;
	extern Residue C;
	extern Residue D;
	extern Residue E;
	extern Residue F;
	extern Residue G;
	extern Residue H;
	extern Residue I;
	extern Residue K;
	extern Residue L;
	extern Residue M;
	extern Residue N;
	extern Residue P;
	extern Residue Q;
	extern Residue R;
	extern Residue S;
	extern Residue T;
	extern Residue V;
	extern Residue W;
	extern Residue Y;

	extern std::map<const char, const Residue*> name2residue;


}


/*

struct NeutralLosses {
	 constexpr Residue H2O  { 71.037114,    2, 0, 1, 0, 0};
	 constexpr Residue NH3  { 103.009185,   3, 1, 0, 0, 0};

	constexpr NeutralLosses() {};
};
*/


#endif //MSACQUISITIONSIMULATOR_RESIDUE_H
