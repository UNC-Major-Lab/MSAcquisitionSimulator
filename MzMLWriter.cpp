//
// Created by Dennis Goldfarb on 9/4/15.
//


#include "MzMLWriter.h"

void MzMLWriter::close_file() {
	out.close();
}


void MzMLWriter::output_file_start() {
	out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
	out << "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\" version=\"1.1.0\" id=\"ms1_and_ms2\">" << std::endl;
	out << "\t<cvList count=\"3\">" << std::endl;
	out << "\t\t<cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" URI=\"http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\" version=\"3.29.0\"/>" << std::endl;
	out << "\t\t<cv id=\"UO\" fullName=\"Unit Ontology\" URI=\"http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo\" version=\"12:10:2011\"/>" << std::endl;
	out << "\t\t<cv id=\"IMS\" fullName=\"Imaging MS Ontology\" URI=\"http://www.maldi-msi.org/download/imzml/imagingMS.obo\" version=\"0.9.1\"/>" << std::endl;
	out << "\t</cvList>" << std::endl;
	out << "\t<softwareList count=\"1\">" << std::endl;
	out << "\t\t<software id=\"JAMSS\" version=\"beta\">" << std::endl;
	out << "\t\t</software>" << std::endl;
	out << "\t</softwareList>" << std::endl;
	out << "\t<instrumentConfigurationList count=\"2\">" << std::endl;
	out << "\t\t<instrumentConfiguration id=\"IC\">" << std::endl;
	out << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000031\" name=\"instrument model\"/>" << std::endl;
	out << "\t\t\t<componentList count=\"0\">" << std::endl;
	out << "\t\t\t</componentList>" << std::endl;
	out << "\t\t</instrumentConfiguration>" << std::endl;
	out << "\t\t<instrumentConfiguration id=\"IC2\">" << std::endl;
	out << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000556\" name=\"LTQ Orbitrap XL\" value=\"\"/>" << std::endl;
	out << "\t\t\t<componentList count=\"3\">" << std::endl;
	out << "\t\t\t\t<source order=\"1\">" << std::endl;
	out << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000073\" name=\"electrospray ionization\" value=\"\"/>" << std::endl;
	out << "\t\t\t\t</source>" << std::endl;
	out << "\t\t\t\t<analyzer order=\"1\">" << std::endl;
	out << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000083\" name=\"radial ejection linear ion trap\" value=\"\"/>" << std::endl;
	out << "\t\t\t\t</analyzer>" << std::endl;
	out << "\t\t\t\t<detector order=\"1\">" << std::endl;
	out << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000253\" name=\"electron multiplier\" value=\"\"/>" << std::endl;
	out << "\t\t\t\t</detector>" << std::endl;
	out << "\t\t\t</componentList>" << std::endl;
	out << "\t\t<softwareRef ref=\"JAMSS\"/>" << std::endl;
	out << "\t\t</instrumentConfiguration>" << std::endl;
	out << "\t</instrumentConfigurationList>" << std::endl;
	out << "\t<dataProcessingList count=\"1\">" << std::endl;
	out << "\t\t<dataProcessing id=\"did_nothing\">" << std::endl;
	out << "\t\t</dataProcessing>" << std::endl;
	out << "\t</dataProcessingList>" << std::endl;
	out << "\t<run id=\"simulated_run\" defaultInstrumentConfigurationRef=\"IC\">" << std::endl;
	//int spectrumCount = totalScans + totalScans * highestNMS2;
	//out << "\t\t<spectrumList count=\"" + spectrumCount + "\" defaultDataProcessingRef=\"did_nothing\">" << std::endl;
	out << "\t\t<spectrumList defaultDataProcessingRef=\"did_nothing\">" << std::endl;
	out.flush();
}

void MzMLWriter::add_to_write_buffer(std::unique_ptr<Scan> s) {
	buffer.push_back(std::move(s));
}

void MzMLWriter::write_buffer() {
	for (auto scan = buffer.begin(); scan != buffer.end(); ++scan) {
		write_scan(scan->get());
	}
	buffer.clear();
}


void MzMLWriter::write_scan(Scan *s) {
	switch (s->scan_type) {
		case Scan::ScanType::MS1:
			write_scan_ms1(static_cast<MS1Scan*>(s));
			break;
		case Scan::ScanType::MS2:
			write_scan_ms2(static_cast<MS2Scan*>(s));
			break;
	}
}


void MzMLWriter::output_file_end() {
	out << "\t\t</spectrumList>" << std::endl;
	out << "\t</run>" << std::endl;
	out << "</mzML>" << std::endl;
	out.flush();
	write_index();
	close_file();
}

void MzMLWriter::write_index() {
	long index_offset = out.tellp(); //TODO TEST THIS
	out << "<indexList count=\"" << offsets.size() << "\">" << std::endl;
	out << "<index name=\"spectrum\">" << std::endl;

	for (int i = 0; i < offsets.size(); i++) {
		int spectrum_index = i+1;
		out << "<offset idRef=\"" << spectrum_index << "\" nativeID=\"" << spectrum_index << "\">" << offsets[i] << "</offset>" << std::endl;
	}
	out << "</index>" << std::endl;
	out << "</indexList>" << std::endl;
	out << "<indexListOffset>" << index_offset << "</indexListOffset>" << std::endl;

}

void MzMLWriter::write_scan_ms1(MS1Scan *s) {
	std::string encodedMZ = compress_peaks(s->peaks, true);
	std::string encodedINT = compress_peaks(s->peaks, false);
		
	offsets.push_back(out.tellp());
	out << "\t\t\t<spectrum index=\"" << s->scan_id << "\" id=\"" << s->scan_id << "\" defaultArrayLength=\"" << s->peaks.size() << "\">" << std::endl;
	out << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\"/>" << std::endl;
	out << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"1\"/>" << std::endl;
	out << "\t\t\t\t<scanList count=\"1\">" << std::endl;
	out << "\t\t\t\t\t<scan>" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"" << s->retention_time << "\" unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"/>" << std::endl;
	out << "\t\t\t\t\t</scan>" << std::endl;
	out << "\t\t\t\t</scanList>" << std::endl;
	out << "\t\t\t\t<binaryDataArrayList count=\"2\">" << std::endl;
	out << "\t\t\t\t\t<binaryDataArray encodedLength=\"" << encodedMZ.length() << "\">" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\"/>" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\"/>" << std::endl;
	out << "\t\t\t\t\t\t<binary>" + encodedMZ + "</binary>" << std::endl;
	out << "\t\t\t\t\t</binaryDataArray>" << std::endl;
	out << "\t\t\t\t\t<binaryDataArray encodedLength=\"" << encodedINT.length() << "\">" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\"/>" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\"/>" << std::endl;
	out << "\t\t\t\t\t\t<binary>" + encodedINT + "</binary>" << std::endl;
	out << "\t\t\t\t\t</binaryDataArray>" << std::endl;
	out << "\t\t\t\t</binaryDataArrayList>" << std::endl;
	out << "\t\t\t</spectrum>" << std::endl;
	out.flush();
}

void MzMLWriter::write_scan_ms2(MS2Scan *s) {
	if (s->peaks.size() == 0) return;

	offsets.push_back(out.tellp());
	std::string encodedMS2MZ = compress_peaks(s->peaks, true);
	std::string encodedMS2INT = compress_peaks(s->peaks, false);
	out << "\t\t\t<spectrum index=\"" << s->scan_id << "\" id=\"" << s->scan_id << "\" defaultArrayLength=\"" << s->peaks.size() << "\">" << std::endl;
	out << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"2\"/>" << std::endl;
	out << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000130\" name=\"positive scan\" value=\"\"/>" << std::endl;
	out << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>" << std::endl;
	out << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\"/>" << std::endl;
	out << "\t\t\t\t<scanList count=\"1\">" << std::endl;
	out << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\" value=\"\"/>" << std::endl;
	out << "\t\t\t\t\t<scan instrumentConfigurationRef=\"IC2\">" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"" << s->retention_time << "\" unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"/>" << std::endl;
	out << "\t\t\t\t\t</scan>" << std::endl;
	out << "\t\t\t\t</scanList>" << std::endl;
	out << "\t\t\t\t<precursorList count=\"1\">" << std::endl;
	out << "\t\t\t\t\t<precursor spectrumRef=\"scan=" << s->parent_scan_id << "\">" << std::endl;
	out << "\t\t\t\t\t\t<selectedIonList count=\"1\">" << std::endl;
	out << "\t\t\t\t\t\t\t<selectedIon>" << std::endl;
	out << "\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000042\" name=\"peak intensity\" value=\"" << s->precursor_peak.intensity << "\"/>" << std::endl;
	out << "\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000041\" name=\"charge state\" value=\"" << "2" <<  "\"/>" << std::endl;
	out << "\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000744\" name=\"selected ion m/z\" value=\"" << s->precursor_peak.mz << "\"/>" << std::endl;
	out << "\t\t\t\t\t\t\t</selectedIon>" << std::endl;
	out << "\t\t\t\t\t\t</selectedIonList>" << std::endl;
	out << "\t\t\t\t\t\t<activation>" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000133\" name=\"collision-induced dissociation\" value=\"\"/>";
	out << "\t\t\t\t\t\t</activation>" << std::endl;
	out << "\t\t\t\t\t</precursor>" << std::endl;
	out << "\t\t\t\t</precursorList>" << std::endl;
	out << "\t\t\t\t<binaryDataArrayList count=\"2\">" << std::endl;
	out << "\t\t\t\t\t<binaryDataArray encodedLength=\"" << encodedMS2MZ.length() << "\">" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\"/>" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\"/>" << std::endl;
	out << "\t\t\t\t\t\t<binary>" + encodedMS2MZ + "</binary>" << std::endl;
	out << "\t\t\t\t\t</binaryDataArray>" << std::endl;
	out << "\t\t\t\t\t<binaryDataArray encodedLength=\"" << encodedMS2INT.length() << "\">" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\"/>" << std::endl;
	out << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\"/>" << std::endl;
	out << "\t\t\t\t\t\t<binary>" + encodedMS2INT + "</binary>" << std::endl;
	out << "\t\t\t\t\t</binaryDataArray>" << std::endl;
	out << "\t\t\t\t</binaryDataArrayList>" << std::endl;
	out << "\t\t\t</spectrum>" << std::endl;
	out.flush();
}

std::string MzMLWriter::compress_peaks(std::vector<BasicPeak> &peaks, bool is_mz) {
	double array_peaks[peaks.size()];

	for (int i = 0; i < peaks.size(); i++) {
		array_peaks[i] = is_mz ? peaks[i].mz : peaks[i].intensity;
	}

	// Compress using zlib
	unsigned char *zCompr=NULL;
	uLong len, zlen;

	len = (uLong) peaks.size()*sizeof(double);
	zlen = compressBound(len);
	zCompr = (unsigned char*)calloc((uInt) zlen, 1);
	compress(zCompr, &zlen, (const Bytef*)&array_peaks[0], len);

	return base64_encode(zCompr, zlen);

}
