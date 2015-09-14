//
// Created by Dennis Goldfarb on 6/29/15.
//

#ifndef MSACQUISITIONSIMULATOR_ELEMENT_H
#define MSACQUISITIONSIMULATOR_ELEMENT_H


#include <iostream>
#include <map>

namespace elements {

	// This order must match the order in libmercury++.h
	enum ELEMENTS {
		H, C, N, O, S, Br, Cl, I, P, F, Si, Fe, B, Se, Co, Mg, Au, As, Sn, Na, K, Te, Zn, Ge, Ca, Sb, Cu, Al, Mn,
		Pt, Gd, Hg, Mo, Sr, Ga, Ni, Pb, Ag, Bi, Tl, Cr, Rb, Zr, Ti, W, Be, V, Cd, Ba, Ta, Li, Cs, Pd, Ce, Ru, La,
		Nd, Re, Hf, Th, He, Ar, Lu, U, Kr, Ir, In, Rh, Ho, Dy, Yb, Eu, Os, Pr, Tb, Er, Xe, Sc, Ne, Sm, Tm, Nb
	};

	const std::map<const std::string, ELEMENTS> name2element{{"H",  H}, {"C",  C}, {"N",  N}, {"O",  O}, {"S",  S}, {"Br",Br}, {"Cl",Cl}, {"I",  I}, {"P",  P},
															 {"F",  F}, {"Si",Si}, {"Fe",Fe}, {"B",  B}, {"Se",Se}, {"Co",Co}, {"Mg",Mg}, {"Au",Au}, {"As",As},
															 {"Sn",Sn}, {"Na",Na}, {"K",  K}, {"Te",Te}, {"Zn",Zn}, {"Ge",Ge}, {"Ca",Ca}, {"Sb",Sb}, {"Cu",Cu},
															 {"Al",Al}, {"Mn",Mn}, {"Pt",Pt}, {"Gd",Gd}, {"Hg",Hg}, {"Mo",Mo}, {"Sr",Sr}, {"Ga",Ga}, {"Ni",Ni},
															 {"Pb",Pb}, {"Ag",Ag}, {"Bi",Bi}, {"Tl",Tl}, {"Cr",Cr}, {"Rb",Rb}, {"Zr",Zr}, {"Ti",Ti}, {"W",  W},
															 {"Be",Be}, {"V",  V}, {"Cd",Cd}, {"Ba",Ba}, {"Ta",Ta}, {"Li",Li}, {"Cs",Cs}, {"Pd",Pd}, {"Ce",Ce},
															 {"Ru",Ru}, {"La",La}, {"Nd",Nd}, {"Re",Re}, {"Hf",Hf}, {"Th",Th}, {"He",He}, {"Ar",Ar}, {"Lu",Lu},
															 {"U",  U}, {"Kr",Kr}, {"Ir",Ir}, {"In",In}, {"Rh",Rh}, {"Ho",Ho}, {"Dy",Dy}, {"Yb",Yb}, {"Eu",Eu},
															 {"Os",Os}, {"Pr",Pr}, {"Tb",Tb}, {"Er",Er}, {"Xe",Xe}, {"Sc",Sc}, {"Ne",Ne}, {"Sm",Sm}, {"Tm",Tm},
															 {"Nb",Nb}};
}
#endif //MSACQUISITIONSIMULATOR_ELEMENT_H
