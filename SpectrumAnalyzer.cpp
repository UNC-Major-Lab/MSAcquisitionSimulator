//
// Created by Dennis Goldfarb on 8/28/15.
//

#include "SpectrumAnalyzer.h"


void SpectrumAnalyzer::deisotope(std::vector<BasicPeak> &peaks) {
	double threshold = 0.01;
	//std::cout << "before: " << peaks.size() << std::endl;
	for (auto mono_peak = peaks.begin(); mono_peak < peaks.end();) {
		for (int z = 1; z <= 5; z++) {
			double step = 1.0/z;
			double m1_mz = mono_peak->mz + step;
			double m2_mz = mono_peak->mz + 2*step;

			auto m1_peak = std::next(mono_peak);
			bool m1_found = false;

			for (; m1_peak != peaks.end() && m1_peak->mz - m1_mz <= threshold; ++m1_peak) {
				if (std::abs(m1_peak->mz - m1_mz) <= threshold && m1_peak->intensity < 1.2*mono_peak->intensity) {
					m1_found = true;
					break;
				}
			}

			if (m1_peak == peaks.end()) continue;
			auto m2_peak = std::next(m1_peak);
			bool m2_found = false;

			for (; m2_peak != peaks.end() && m2_peak->mz - m2_mz <= threshold; ++m2_peak) {
				if (std::abs(m2_peak->mz - m2_mz) <= threshold && m2_peak->intensity < m1_peak->intensity) {
					m2_found = true;
					break;
				}
			}

			if (m1_found && m2_found) {
				double next_mz = m2_mz+step;
				double prev_int = m2_peak->intensity;
				if (m2_peak != peaks.end()) {
					for (auto next_peak = std::next(m2_peak); next_peak != peaks.end() && next_peak->mz - next_mz <= threshold;) {
						if (std::abs(next_peak->mz - next_mz) <= threshold && next_peak->intensity < prev_int) {
							next_mz += step;
							prev_int = next_peak->intensity;
							next_peak = peaks.erase(next_peak);
						} else {
							++next_peak;
						}
					}
				}
				peaks.erase(m1_peak);
				peaks.erase(m2_peak);
				break;
			}
		}

		++mono_peak;
	}
	//std::cout << "after: " << peaks.size() << std::endl;
	return;
}



