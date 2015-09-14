//
// Created by Dennis Goldfarb on 9/4/15.
//

#ifndef MSACQUISITIONSIMULATOR_MZXMLWRITER_H
#define MSACQUISITIONSIMULATOR_MZXMLWRITER_H


#include <array>
#include <fstream>
#include <vector>
#include <memory>
#include <zlib.h>
#include "Scan.h"
#include "MS1Scan.h"
#include "MS2Scan.h"
#include "base64.h"

class MzMLWriter {
private:
	std::ofstream out;
	std::vector<long> offsets;


	void write_scan(Scan* s);
	void write_scan_ms1(MS1Scan* s);
	void write_scan_ms2(MS2Scan* s);
	void output_file_start();

	void write_index();
	std::string compress_peaks(std::vector<BasicPeak> &peaks, bool is_mz);

public:
	MzMLWriter(std::string output_path) : out(output_path.c_str()) {
		output_file_start();
	};

	void close_file();
	void add_to_scan_buffer(std::unique_ptr<Scan> s);
	void write_buffer();
	void output_file_end();

	std::vector<std::unique_ptr<Scan>> buffer;
};


#endif //MSACQUISITIONSIMULATOR_MZXMLWRITER_H
