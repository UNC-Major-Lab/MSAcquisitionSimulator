//
// Created by Dennis Goldfarb on 9/9/15.
//

#include "FidoWriter.h"

void FidoWriter::write_peptide(double probablity, std::string peptide, std::vector<Protein *> proteins) {
	out << "e " << peptide << std::endl;
	for (Protein* p : proteins) {
		out << "r " << p->name << std::endl;
	}
	out << "p " << probablity << std::endl;
}

void FidoWriter::close_file() {
	out.close();
}
