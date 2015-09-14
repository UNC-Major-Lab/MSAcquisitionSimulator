//
// Created by Dennis Goldfarb on 6/29/15.
//

#include "Residue.h"


namespace residues {

	Residue A = Residue {71.037114,    "H5,N1,C3,O1"};
	Residue C = Residue {103.009185,   "H5,N1,C3,O1,S1"};
	Residue D = Residue {115.026943,   "H5,N1,C4,O3"};
	Residue E = Residue {129.042593,   "H7,N1,C5,O3"};
	Residue F = Residue {147.068414,   "H9,N1,C9,O1"};
	Residue G = Residue {57.021464,    "H3,N1,C2,O1"};
	Residue H = Residue {137.058912,   "H7,N3,C6,O1"};
	Residue I = Residue {113.084064,   "H11,N1,C6,O1"};
	Residue K = Residue {128.094963,   "H12,N2,C6,O1"};
	Residue L = Residue {113.084064,   "H11,N1,C6,O1"};
	Residue M = Residue {131.040485,   "H9,N1,C5,O1,S1"};
	Residue N = Residue {114.042927,   "H6,N2,C4,O2"};
	Residue P = Residue {97.052764,    "H7,N1,C5,O1"};
	Residue Q = Residue {128.058578,   "H8,N2,C5,O2"};
	Residue R = Residue {156.101111,   "H12,N4,C6,O1"};
	Residue S = Residue {87.032028,    "H5,N1,C3,O2"};
	Residue T = Residue {101.047679,   "H7,N1,C4,O2"};
	Residue V = Residue {99.068414,    "H9,N1,C5,O1"};
	Residue W = Residue {186.079313,   "H10,N2,C11,O1"};
	Residue Y = Residue {163.06332,    "H9,N1,C9,O2"};

	std::map<const char, const Residue*> name2residue = std::map<const char, const Residue*>{{'A', &A}, {'C', &C}, {'D', &D}, {'E', &E},
																							 {'F', &F}, {'G', &G}, {'H', &H}, {'I', &I},
																							 {'K', &K}, {'L', &L}, {'M', &M}, {'N', &N},
																							 {'P', &P}, {'Q', &Q}, {'R', &R}, {'S', &S},
																							 {'T', &T}, {'V', &V}, {'W', &W}, {'Y', &Y}};


}