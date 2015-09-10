//
// Created by Dennis Goldfarb on 8/12/15.
//

#include "GroundTruthDatabase.h"
#include "Globals.h"

void GroundTruthDatabase::initialize_database() {
	if (create_new) {
		drop_database();
		create_database();
		create_tables();
		create_indexes();
	} else {
		open_database();
	}
	prepare_stmts();
}

void GroundTruthDatabase::open_database() {
	sqlite3_backup *pBackup;
	sqlite3 *pFrom;

	int rc;
	rc = sqlite3_open(path.c_str(), &pFrom);
	rc = sqlite3_open(":memory:", &db);
	pBackup = sqlite3_backup_init(db, "main", pFrom, "main");
	if( pBackup ){
		sqlite3_backup_step(pBackup, -1);
		sqlite3_backup_finish(pBackup);
	}

	sqlite3_close(pFrom);

	if (rc) {
		std::cerr << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
		sqlite3_close(db);
	}
}

int GroundTruthDatabase::create_database() {
	int rc;
	rc = sqlite3_open(path.c_str(), &db);
	if (rc) {
		std::cerr << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
		sqlite3_close(db);
		return 1;
	}
	sqlite3_exec(db, "PRAGMA synchronous = OFF", NULL, NULL, NULL);
	sqlite3_exec(db, "PRAGMA journal_mode = MEMORY", NULL, NULL, NULL);
	return 0;
}

int GroundTruthDatabase::close_database() {
	finalize_stmts();
	int rc = sqlite3_close(db);
	return rc;
}

int GroundTruthDatabase::create_tables() {
	execute_sql("CREATE TABLE proteins ("
			"abundance REAL NOT NULL, "
			"sequence TEXT NOT NULL, "
			"protein_name TEXT NOT NULL);");

	execute_sql("CREATE TABLE peptides ("
			"protein_id INT NOT NULL, "
			"start_index INT NOT NULL, "
			"end_index INT NOT NULL);");

	execute_sql("CREATE TABLE mods ("
			"name TEXT NOT NULL, "
			"abbreviation TEXT NOT NULL);");

	execute_sql("CREATE TABLE peptide_mods ("
			"aa_index INT NOT NULL, "
			"mod_id INT NOT NULL, "
			"peptide_id INT NOT NULL,"
			"PRIMARY KEY (aa_index, mod_id, peptide_id));");

	execute_sql("CREATE TABLE ions ("
			"peptide_id INT NOT NULL, "
			"neutrons INT NOT NULL, "
			"charge INT NOT NULL, "
			"rt REAL NOT NULL, "
			"rt_start REAL NOT NULL, "
			"rt_end REAL NOT NULL, "
			"intensity REAL NOT NULL, "
			"mz REAL NOT NULL);");

	return 0;
}

int GroundTruthDatabase::execute_sql(std::string sql) {
	int rc;
	sqlite3_stmt *stmt;
	rc = sqlite3_prepare(db, sql.c_str(), -1, &stmt, 0);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
	return rc;
}

void GroundTruthDatabase::insert_ions(const double &abundance, const int charge, const double rt, const long peptide_id,
									 const std::vector<double> &isotope_mz,
									 const std::vector<double> &isotope_abundance, ElutionShapeSimulator &elution_shape_simulator) {
	for (int i = 0; i < isotope_mz.size(); ++i) {
		double mz = isotope_mz[i];
		if (mz < MIN_MZ || mz > MAX_MZ) continue;
		double iso_abundance = isotope_abundance[i] * abundance;
		if (iso_abundance < PRUNE_THRESHOLD) break;
		iso_abundance /= elution_shape_simulator.normalization_factor;

		double rt_start = elution_shape_simulator.get_min_rt(rt,iso_abundance);
		double rt_end = elution_shape_simulator.get_max_rt(rt,iso_abundance);

		if (std::isnan(rt_start) || std::isnan(rt_end)) continue;

		sqlite3_bind_int64(insert_ion_stmt, 1, peptide_id);
		sqlite3_bind_int(insert_ion_stmt, 2, i);
		sqlite3_bind_int(insert_ion_stmt, 3, charge);
		sqlite3_bind_double(insert_ion_stmt, 4, rt);
		sqlite3_bind_double(insert_ion_stmt, 5, rt_start);
		sqlite3_bind_double(insert_ion_stmt, 6, rt_end);
		sqlite3_bind_double(insert_ion_stmt, 7, iso_abundance);
		sqlite3_bind_double(insert_ion_stmt, 8, mz);

		if(sqlite3_step(insert_ion_stmt) != SQLITE_DONE) {
			std::cout << sqlite3_errmsg(db) << std::endl;
		}
		sqlite3_reset(insert_ion_stmt);
	}
}

long GroundTruthDatabase::insert_peptide(const Peptide &p, const long protein_id, const std::map<PTM*,long> &ptm2id) {
	sqlite3_bind_int64(insert_peptide_stmt, 1, protein_id);
	sqlite3_bind_int(insert_peptide_stmt, 2, p.start);
	sqlite3_bind_int(insert_peptide_stmt, 3, p.end);

	if(sqlite3_step(insert_peptide_stmt) != SQLITE_DONE) {
		std::cout << sqlite3_errmsg(db) << std::endl;
	}

	long peptide_id = sqlite3_last_insert_rowid(db);
	sqlite3_reset(insert_peptide_stmt);

	for (std::pair<int,PTM*> pair : p.index2mod_for_pep) {
		sqlite3_bind_int(insert_peptide_modification_stmt, 1, pair.first);
		sqlite3_bind_int64(insert_peptide_modification_stmt, 2, ptm2id.at(pair.second));
		sqlite3_bind_int64(insert_peptide_modification_stmt, 3, peptide_id);

		int rc = sqlite3_step(insert_peptide_modification_stmt);
		if(rc != SQLITE_DONE) {
			std::cout << "insert peptide modification: " << rc << " " << sqlite3_errmsg(db) << std::endl;
		}

		sqlite3_reset(insert_peptide_modification_stmt);
	}

	return peptide_id;
}

long GroundTruthDatabase::insert_protein(const Protein &p) {
	sqlite3_bind_double(insert_protein_stmt, 1, p.abundance);
	sqlite3_bind_text(insert_protein_stmt, 2, p.sequence.c_str(), (int) p.sequence.length()+1, SQLITE_TRANSIENT);
	sqlite3_bind_text(insert_protein_stmt, 3, p.name.c_str(), (int) p.name.length()+1, SQLITE_TRANSIENT);

	if(sqlite3_step(insert_protein_stmt) != SQLITE_DONE) {
		std::cout << sqlite3_errmsg(db) << std::endl;
	}

	long protein_id = sqlite3_last_insert_rowid(db);
	sqlite3_reset(insert_protein_stmt);

	return protein_id;
}

long GroundTruthDatabase::insert_modification(const PTM &ptm) {
	sqlite3_bind_text(insert_modification_stmt, 1, ptm.name.c_str(), (int) ptm.name.length()+1, SQLITE_TRANSIENT);
	sqlite3_bind_text(insert_modification_stmt, 2, ptm.abbreviation.c_str(), (int) ptm.abbreviation.length()+1, SQLITE_TRANSIENT);

	if(sqlite3_step(insert_modification_stmt) != SQLITE_DONE) {
		std::cout << sqlite3_errmsg(db) << std::endl;
	}

	long mod_id = sqlite3_last_insert_rowid(db);
	sqlite3_reset(insert_modification_stmt);

	return mod_id;
}

void GroundTruthDatabase::drop_database() {
	std::ifstream infile(path.c_str());
	if (infile) {
		if (remove(path.c_str()) != 0)
			std::cerr << "Error deleting file" << std::endl;
	}
}

void GroundTruthDatabase::prepare_stmts() {
	if (create_new) {
		sqlite3_prepare(db, "insert into mods (name, abbreviation) values (?,?)", -1, &insert_modification_stmt, 0);
		sqlite3_prepare(db, "insert into proteins (abundance, sequence, protein_name) values (?,?,?)", -1,&insert_protein_stmt, 0);
		sqlite3_prepare(db, "insert into peptides (protein_id, start_index, end_index) values (?,?,?)", -1, &insert_peptide_stmt, 0);
		sqlite3_prepare(db, "insert into peptide_mods (aa_index, mod_id, peptide_id) values (?,?,?)", -1, &insert_peptide_modification_stmt, 0);
		sqlite3_prepare(db, "insert into ions (peptide_id, neutrons, charge, rt, rt_start, rt_end, intensity, mz) values (?,?,?,?,?,?,?,?)", -1, &insert_ion_stmt, 0);
	} else {
		sqlite3_prepare(db, "select neutrons, charge, rt, intensity, mz from ions where rt_start <= ? and rt_end >= ? and mz >= ? and mz <= ?", -1, &select_ions_stmt, 0);
		sqlite3_prepare(db, "select i.ROWID, neutrons, charge, rt, intensity, mz, prot.sequence, pep.start_index, pep.end_index, "
				"pmods.aa_index, m.abbreviation "
				"from ions i inner join peptides pep on pep.ROWID = i.peptide_id "
				"left join peptide_mods pmods on pmods.peptide_id = pep.ROWID "
				"left join mods m on m.ROWID = pmods.mod_id "
				"inner join proteins prot on prot.ROWID = pep.protein_id "
				"where rt_start <= ? and rt_end >= ? and mz >= ? and mz <= ? "
				"order by i.ROWID", -1, &select_peptides_stmt, 0);
	}
}

void GroundTruthDatabase::finalize_stmts() {
	if (create_new) {
		sqlite3_finalize(insert_modification_stmt);
		sqlite3_finalize(insert_protein_stmt);
		sqlite3_finalize(insert_peptide_stmt);
		sqlite3_finalize(insert_peptide_modification_stmt);
		sqlite3_finalize(insert_ion_stmt);
	} else {
		sqlite3_finalize(select_ions_stmt);
		sqlite3_finalize(select_peptides_stmt);
	}
}

void GroundTruthDatabase::create_indexes() {
	sqlite3_exec(db, "CREATE INDEX retention_time ON ions (rt_start, rt_end, mz)", 0, 0, 0);
}


std::vector<Centroid*> GroundTruthDatabase::get_ions_at_rt(double min_mz, double max_mz, double time) {
	sqlite3_bind_double(select_ions_stmt, 1, time);
	sqlite3_bind_double(select_ions_stmt, 2, time);
	sqlite3_bind_double(select_ions_stmt, 3, min_mz);
	sqlite3_bind_double(select_ions_stmt, 4, max_mz);

	std::vector<Centroid*> ions;

	int rc;
	while ((rc = sqlite3_step(select_ions_stmt)) == SQLITE_ROW) {
		int num_neutrons = sqlite3_column_int(select_ions_stmt, 0);
		int charge = sqlite3_column_int(select_ions_stmt, 1);
		double rt_center = sqlite3_column_double(select_ions_stmt, 2);
		double intensity = sqlite3_column_double(select_ions_stmt, 3);
		double mz = sqlite3_column_double(select_ions_stmt, 4);
		double mass = mz * charge;
		ions.push_back(new Centroid(mz, mass, charge, intensity, num_neutrons, rt_center));
	}

	sqlite3_reset(select_ions_stmt);

	return ions;
}

std::vector<Centroid*> GroundTruthDatabase::get_peptides_at_rt(double min_mz, double max_mz, double time) {
	sqlite3_bind_double(select_peptides_stmt, 1, time);
	sqlite3_bind_double(select_peptides_stmt, 2, time);
	sqlite3_bind_double(select_peptides_stmt, 3, min_mz);
	sqlite3_bind_double(select_peptides_stmt, 4, max_mz);

	std::vector<Centroid*> ions;

	int rc;

	int prev_id = -1;
	int num_neutrons, charge, pep_start, pep_end, mod_index;
	double rt_center, intensity, mz, mass;
	const unsigned char* prot_seq;
	const unsigned char* abbr;
	std::string sequence;
	std::map<int,std::string> index2mod;

	while ((rc = sqlite3_step(select_peptides_stmt)) == SQLITE_ROW) {
		int ion_id = sqlite3_column_int(select_peptides_stmt, 0);

		if (ion_id == prev_id) { // add PTM
			mod_index = sqlite3_column_int(select_peptides_stmt, 9);
			abbr = sqlite3_column_text(select_peptides_stmt, 10);

			if (abbr) index2mod[mod_index] = std::string(reinterpret_cast<const char*>(abbr));

		} else {
			if (prev_id >= 0) { // append MS2Centroid
				std::string modified_sequence;

				for (int i = pep_start; i < sequence.size()+pep_start; ++i) {
					if (index2mod.find(i) != index2mod.end()) {
						modified_sequence += index2mod.at(i);
					}
					modified_sequence.push_back(sequence[i-pep_start]);
				}

				ions.push_back(new MS2Centroid(mz, mass, charge, intensity, num_neutrons, rt_center, modified_sequence));
				index2mod.clear();
			}
			num_neutrons = sqlite3_column_int(select_peptides_stmt, 1);
			charge = sqlite3_column_int(select_peptides_stmt, 2);
			rt_center = sqlite3_column_double(select_peptides_stmt, 3);
			intensity = sqlite3_column_double(select_peptides_stmt, 4);
			mz = sqlite3_column_double(select_peptides_stmt, 5);
			prot_seq = sqlite3_column_text(select_peptides_stmt, 6);
			pep_start = sqlite3_column_int(select_peptides_stmt, 7);
			pep_end = sqlite3_column_int(select_peptides_stmt, 8);
			mod_index = sqlite3_column_int(select_peptides_stmt, 9);
			abbr = sqlite3_column_text(select_peptides_stmt, 10);
			mass = mz * charge;

			if (abbr) index2mod[mod_index] = std::string(reinterpret_cast<const char*>(abbr));

			sequence = std::string(prot_seq+pep_start, prot_seq+pep_end+1);

			prev_id = ion_id;
		}

	}

	sqlite3_reset(select_peptides_stmt);

	return ions;
}
