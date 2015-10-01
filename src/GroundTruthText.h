//
// Created by Dennis Goldfarb on 9/10/15.
//

#ifndef MSACQUISITIONSIMULATOR_GROUNDTRUTHTEXT_H
#define MSACQUISITIONSIMULATOR_GROUNDTRUTHTEXT_H

#include <fstream>
#include <vector>
#include <string>
#include <list>
#include "IonDAO.h"
#include "Peptide.h"
#include "ElutionShapeSimulator.h"

class GroundTruthText {
private:
	void close_file();
	void insert_ion(IonDAO &ion);
	std::list<IonDAO*> current_ions;

public:

	GroundTruthText(std::string path, bool create_new);

	~GroundTruthText() {
		close_file();
	}

	void insert_ions(double abundance, int charge, double rt, const Peptide* peptide,
					 std::vector<double> &isotope_mz, std::vector<double> &isotope_abundance,
					 ElutionShapeSimulator &elution_shape_simulato, int &num_processed_ionsr);

	std::ofstream *out;
	std::ifstream *in;
	std::string path;
	bool create_new;
	void write_sorted_file();

	std::vector<IonDAO*> get_ions_at_rt(double min_mz, double max_mz, double time);
};


#endif //MSACQUISITIONSIMULATOR_GROUNDTRUTHTEXT_H
