//
// Created by Dennis Goldfarb on 6/27/15.
//

#include "FastaParser.h"
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

std::vector<Protein> parse_FASTA(const char * path_FASTA) {
	std::vector<Protein> proteins = std::vector<Protein>();
	gzFile fp;
	fp = gzopen(path_FASTA, "r"); // STEP 2: open the file handler
	if (!fp) {
		std::cerr << "Error opening FASTA file: " << path_FASTA << std::endl;
		throw std::exception();
	}

	kseq_t *seq;
	int l;

	seq = kseq_init(fp); // STEP 3: initialize seq
	while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
		proteins.push_back(Protein(seq->name.s, seq->seq.s, seq->comment.s, parse_abundance(seq->comment.s)));
	}
	kseq_destroy(seq);
	gzclose(fp);

	return proteins;
}

double parse_abundance(std::string description) {
	std::size_t found = description.find_last_of("#");
	if (found != std::string::npos) {
		return stod(description.substr(found+1));
	}
	return 0;
}