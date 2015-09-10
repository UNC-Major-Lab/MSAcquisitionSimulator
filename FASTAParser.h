//
// Created by Dennis Goldfarb on 6/27/15.
//

#ifndef MSACQUISITIONSIMULATOR_FASTAPARSER_H
#define MSACQUISITIONSIMULATOR_FASTAPARSER_H


#include <iostream>
#include <vector>
#include "Protein.h"

double parse_abundance(std::string description);
std::vector<Protein> parse_FASTA(const char * path_FASTA);


#endif //MSACQUISITIONSIMULATOR_FASTAPARSER_H
