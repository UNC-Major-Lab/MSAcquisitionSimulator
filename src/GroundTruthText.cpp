//
// Created by Dennis Goldfarb on 9/10/15.
//

#include "GroundTruthText.h"

GroundTruthText::GroundTruthText(std::string path, bool create_new) : path(path), create_new(create_new) {
	if (create_new) {
		out = new std::ofstream(path);
	} else {
		in = new std::ifstream(path);
	}
}

void GroundTruthText::close_file() {
	if (create_new) {
		out->close();
		delete out;
	} else {
		in->close();
		delete in;
	}
}

int GroundTruthText::insert_ions(double abundance, int charge, double rt, const Peptide *peptide,
								  std::vector<double> &isotope_mz, std::vector<double> &isotope_abundance,
								  ElutionShapeSimulator &elution_shape_simulator) {

	int num_processed_ions = 0;
	for (int i = 0; i < isotope_mz.size(); ++i) {
		double mz = isotope_mz[i];
		if (mz < MIN_MZ || mz > MAX_MZ) continue;
		double iso_abundance = isotope_abundance[i] * abundance;
		if (iso_abundance < PRUNE_THRESHOLD) break;
		iso_abundance /= elution_shape_simulator.normalization_factor;

		double rt_start = elution_shape_simulator.get_min_rt(rt,iso_abundance);
		double rt_end = elution_shape_simulator.get_max_rt(rt,iso_abundance);

		if (std::isnan(rt_start) || std::isnan(rt_end)) continue;

		ions.push_back(IonTinyDAO(i, charge, iso_abundance, rt, rt_start, rt_end, mz, peptide));
		num_processed_ions++;
	}
	return num_processed_ions;
}

void GroundTruthText::write_sorted_file() {
	std::sort(ions.begin(), ions.end(), IonTinyDAO::less_rt_start);

	for (IonTinyDAO &ion : ions) insert_ion(ion);
}

void GroundTruthText::insert_ion(IonTinyDAO &ion) {
	*out << ion.abundance << "\t" << ion.neutrons << "\t" << ion.mz << "\t" << ion.charge << "\t" << ion.rt << "\t"
	<< ion.rt_start << "\t" << ion.rt_end << "\t" << ion.peptide->get_modified_sequence() << "\t" << ion.peptide->start
	<< "\t" << ion.peptide->end << "\t" << std::endl;
}

std::vector<IonDAO *> GroundTruthText::get_ions_at_rt(double min_mz, double max_mz, double time) {
	std::vector<IonDAO *> ions;

	// remove old ions
	for (auto ion = current_ions.begin(); ion != current_ions.end();) {
		if ((*ion)->rt_end < time) {
			delete *ion;
			ion = current_ions.erase(ion);
		} else {
			++ion;
		}
	}

	// get new ions
	double abundance, rt, rt_start, rt_end, mz;
	int neutrons, charge, peptide_start, peptide_end;
	std::string modified_sequence;
	while (!(*in).eof()) {
		*in >> abundance >> neutrons >> mz >> charge >> rt >> rt_start >> rt_end >> modified_sequence >> peptide_start >>
		peptide_end;
		current_ions.push_back(new IonDAO(neutrons, charge, peptide_start, peptide_end, abundance, rt, rt_start, rt_end, mz, modified_sequence));
		if (rt_start > time) break;
	}

	// copy ions over
	for (auto ion = current_ions.begin(); ion != current_ions.end(); ++ion) {
		if ((*ion)->mz >= min_mz && (*ion)->mz <= max_mz) {
			ions.push_back(*ion);
		}
	}

	return ions;
}
