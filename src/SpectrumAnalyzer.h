//
// Created by Dennis Goldfarb on 8/28/15.
//

#ifndef MSACQUISITIONSIMULATOR_SPECTRUMANALYZER_H
#define MSACQUISITIONSIMULATOR_SPECTRUMANALYZER_H

#include <list>
#include <vector>
#include <string>
#include <cmath>
#include "Peak.h"

class SpectrumAnalyzer {
private:


public:
	SpectrumAnalyzer() {};

	void deisotope(std::vector<BasicPeak> &peaks);
};


#endif //MSACQUISITIONSIMULATOR_SPECTRUMANALYZER_H
