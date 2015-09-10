//
// Created by Dennis Goldfarb on 9/9/15.
//

#ifndef MSACQUISITIONSIMULATOR_FIDOWRITER_H
#define MSACQUISITIONSIMULATOR_FIDOWRITER_H

#include <vector>
#include <string>
#include <fstream>
#include "Protein.h"

class FidoWriter {
private:
	std::ofstream out;
public:
	FidoWriter(std::string output_path) : out(output_path.c_str()) {}

	~FidoWriter() {
		close_file();
	}

	void write_peptide(double probablity, std::string peptide, std::vector<Protein*> proteins);
	void close_file();
};


#endif //MSACQUISITIONSIMULATOR_FIDOWRITER_H
