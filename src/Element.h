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
}
#endif //MSACQUISITIONSIMULATOR_ELEMENT_H
