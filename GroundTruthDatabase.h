//
// Created by Dennis Goldfarb on 8/12/15.
//

#ifndef MSACQUISITIONSIMULATOR_GROUNDTRUTHDATABASE_H
#define MSACQUISITIONSIMULATOR_GROUNDTRUTHDATABASE_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <sqlite3.h>
#include "Peptide.h"
#include "PTM.h"
#include "ElutionShapeSimulator.h"
#include "Centroid.h"
#include "MS2Peptide.h"
#include "MS2Centroid.h"

class GroundTruthDatabase {

private:
	int execute_sql(std::string sql);
	int create_tables();
	int create_database();
	void drop_database();
	void prepare_stmts();
	void finalize_stmts();
	void create_indexes();

	void open_database();

	sqlite3_stmt *insert_ion_stmt;
	sqlite3_stmt *insert_protein_stmt;
	sqlite3_stmt *insert_peptide_stmt;
	sqlite3_stmt *insert_modification_stmt;
	sqlite3_stmt *insert_peptide_modification_stmt;

	sqlite3_stmt *select_ions_stmt;
	sqlite3_stmt *select_peptides_stmt;

public:

	GroundTruthDatabase(std::string path, bool create_new) : path(path), create_new(create_new) {
		initialize_database();
	}

	~GroundTruthDatabase() {
		close_database();
	}

	sqlite3 *db;
	std::string path;
	bool create_new;

	void initialize_database();
	int close_database();
	void insert_ions(const double &abundance, const int charge, const double rt, const long peptide_id,
					const std::vector<double> &isotope_mz,
					const std::vector<double> &isotope_abundance, ElutionShapeSimulator &elution_shape_simulator);
	long insert_peptide(const Peptide &p, const long protein_id, const std::map<PTM*,long> &ptm2id);
	long insert_protein(const Protein &p);
	long insert_modification(const PTM &ptm);

	std::vector<Centroid*> get_ions_at_rt(double min_mz, double max_mz, double time);
	std::vector<Centroid*> get_peptides_at_rt(double min_mz, double max_mz, double time);
};


#endif //MSACQUISITIONSIMULATOR_GROUNDTRUTHDATABASE_H
