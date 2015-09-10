#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "biolccc.h"

namespace BioLCCC {

// Auxiliary functions that shouldn't be exposed to a user.
	//namespace {
		double calculateKd(const std::vector<ChemicalGroup> &parsedSequence,
						   const double secondSolventConcentration,
						   const ChemicalBasis &chemBasis,
						   const double columnPoreSize,
						   const double columnRelativeStrength,
						   const double temperature) throw(BioLCCCException) {
			// assymetricCalculations shows whether the Kd for the reversed molecule
			// will differ. It happens when a molecule cannot be divided into an integer
			// number of Kuhn segments.
			bool assymetricCalculations =
					(fmod(chemBasis.monomerLength() * parsedSequence.size(),
						  chemBasis.kuhnLength()) != 0);
			// Choosing the appropriate polymerModel.
			if (chemBasis.polymerModel() == CHAIN) {
				double Kd = calculateKdChain(parsedSequence,
											 secondSolventConcentration,
											 chemBasis, columnPoreSize,
											 columnRelativeStrength, temperature);
				if (assymetricCalculations) {
					std::vector<ChemicalGroup> revParsedSequence =
							parsedSequence;
					std::reverse(revParsedSequence.begin(),
								 revParsedSequence.end());
					Kd = (Kd + calculateKdChain(revParsedSequence,
												secondSolventConcentration,
												chemBasis,
												columnPoreSize,
												columnRelativeStrength,
												temperature)) / 2.0;
				}
				return Kd;
			}
			else if (chemBasis.polymerModel() == ROD) {
				double Kd = calculateKdRod(parsedSequence,
										   secondSolventConcentration,
										   chemBasis, columnPoreSize,
										   columnRelativeStrength, temperature);
				if (assymetricCalculations) {
					std::vector<ChemicalGroup> revParsedSequence =
							parsedSequence;
					std::reverse(revParsedSequence.begin(),
								 revParsedSequence.end());
					Kd = (Kd + calculateKdRod(revParsedSequence,
											  secondSolventConcentration,
											  chemBasis,
											  columnPoreSize,
											  columnRelativeStrength,
											  temperature)) / 2.0;
				}
				return Kd;
			}
			else {
				throw BioLCCCException("Model error.");
			}
		}

		class KdCalculator {
		public:
			KdCalculator(const std::vector<ChemicalGroup> &parsedSequence,
						 const ChemicalBasis &chemBasis,
						 const double columnPoreSize,
						 const double columnRelativeStrength,
						 const double temperature,
						 const int numInterpolationPoints) :
					mParsedSequence(parsedSequence),
					mChemicalBasis(chemBasis),
					mColumnPoreSize(columnPoreSize),
					mColumnRelativeStrength(columnRelativeStrength),
					mTemperature(temperature),
					mN(numInterpolationPoints) {
				if (mN > 0) {
					mSecondSolventConcentrations = new double[mN];
					mLogKds = new double[mN];
					// The number of extra points in the terminal segments.
					// This points significantly increase the accuracy of spline
					// interpolation.
					int NETP = 1;
					for (int i = 0; i < mN; i++) {
						if (i <= NETP) {
							mSecondSolventConcentrations[i] =
									i * 100.0 / (mN - 2.0 * NETP - 1.0) / (NETP + 1.0);
						}
						else if (i > (mN - NETP - 2)) {
							mSecondSolventConcentrations[i] =
									((mN - 2.0 * NETP - 2.0)
									 + (i - mN + NETP + 2.0) / (NETP + 1.0))
									* 100.0 / (mN - 2.0 * NETP - 1.0);
						}
						else {
							mSecondSolventConcentrations[i] =
									(i - NETP) * 100.0 / (mN - 2.0 * NETP - 1.0);
						}
						mLogKds[i] = log(calculateKd(mParsedSequence,
													 mSecondSolventConcentrations[i],
													 mChemicalBasis,
													 mColumnPoreSize,
													 mColumnRelativeStrength,
													 mTemperature));
					}

					mSecondDers = new double[mN];
					fitSpline(mSecondSolventConcentrations, mLogKds,
							  mN, mSecondDers);
				}
			}

			~KdCalculator() {
				if (mN > 0) {
					delete[] mSecondSolventConcentrations;
					delete[] mLogKds;
					delete[] mSecondDers;
				}
			}

			double operator()(double secondSolventConcentration)
			throw(BioLCCCException) {
				if (mN == 0) {
					return calculateKd(mParsedSequence,
									   secondSolventConcentration,
									   mChemicalBasis,
									   mColumnPoreSize,
									   mColumnRelativeStrength,
									   mTemperature);
				}
				else {
					return exp(calculateSpline(mSecondSolventConcentrations, mLogKds,
											   mSecondDers, mN,
											   secondSolventConcentration));
				}
			}

		private:
			const std::vector<ChemicalGroup> &mParsedSequence;
			const ChemicalBasis &mChemicalBasis;
			const double mColumnPoreSize;
			const double mColumnRelativeStrength;
			const double mTemperature;
			const int mN;
			double *mSecondSolventConcentrations;
			double *mLogKds;
			double *mSecondDers;
		};

		double calculateRT(const std::vector<ChemicalGroup> &parsedSequence,
						   const ChemicalBasis &chemBasis,
						   const ChromoConditions &conditions,
						   const int numInterpolationPoints,
						   const bool continueGradient,
						   const bool backwardCompatibility
		) throw(BioLCCCException) {

			if (numInterpolationPoints < 0) {
				throw BioLCCCException(
						"The number of interpolation points must be non-negative.");
			}

			KdCalculator kdCalculator(parsedSequence, chemBasis,
									  conditions.columnPoreSize(),
									  conditions.columnRelativeStrength(),
									  conditions.temperature(),
									  numInterpolationPoints);

			double RT = 0.0;
			// Use simplified expression for isocratic elution.
			if (conditions.SSConcentrations().size() == 1) {
				RT = kdCalculator(conditions.SSConcentrations()[0])
					 * conditions.columnPoreVolume() / conditions.flowRate();
			}
			else {
				// The part of a column passed by molecules. When it exceeds 1.0,
				// the analyte elutes from the column.
				double S = 0.0;
				double dS = 0.0;
				int j = 0;
				double currentSSConcentration = 0.0;
				while (S < 1.0) {
					j++;
					if (j < conditions.SSConcentrations().size()) {
						currentSSConcentration = conditions.SSConcentrations()[j];
					}
					else {
						// If continue gradient then use the slope of the last section.
						if (continueGradient) {
							currentSSConcentration +=
									conditions.SSConcentrations().back()
									- *(conditions.SSConcentrations().end() - 2);
						}
						else {
							break;
						}
					}
					dS = conditions.dV()
						 / kdCalculator(currentSSConcentration)
						 / conditions.columnPoreVolume();
					S += dS;
				}

				RT = j * conditions.dV() / conditions.flowRate();
				// Correction for the discreteness of integration.
				if ((!backwardCompatibility) && (S > 1.0)) {
					RT -= (S - 1.0) / dS * conditions.dV() / conditions.flowRate();
				}
			}

			RT += conditions.delayTime();
			// Correction for V0.
			RT += conditions.columnInterstitialVolume() / conditions.flowRate();
			return RT;
		}
	//}

	double calculateRT(const std::string &sequence,
					   const ChemicalBasis &chemBasis,
					   const ChromoConditions &conditions,
					   const int numInterpolationPoints,
					   const bool continueGradient,
					   const bool backwardCompatibility)
	throw(BioLCCCException) {
		std::vector<ChemicalGroup> parsedSequence =
				parseSequence(sequence, chemBasis);
		return calculateRT(parsedSequence,
						   chemBasis,
						   conditions,
						   numInterpolationPoints,
						   continueGradient,
						   backwardCompatibility);
	}

	double calculateKd(const std::string &sequence, const double secondSolventConcentration,
					   const ChemicalBasis &chemBasis,
					   const double columnPoreSize, const double columnRelativeStrength,
					   const double temperature) throw(BioLCCCException) {

		return calculateKd(parseSequence(sequence, chemBasis), secondSolventConcentration, chemBasis,
						   columnPoreSize, columnRelativeStrength, temperature);
	}

	double calculateAverageMass(const std::string &sequence, const ChemicalBasis &chemBasis) throw(BioLCCCException) {
		std::vector<ChemicalGroup> parsedSequence = parseSequence(sequence, chemBasis);

		double peptideAverageMass = 0;

		for (std::vector<ChemicalGroup>::const_iterator i = parsedSequence.begin(); i < parsedSequence.end(); i++) {
			peptideAverageMass += i->averageMass();
		}

		return peptideAverageMass;
	}

	double calculateMonoisotopicMass(const std::string &sequence,
									 const ChemicalBasis &chemBasis) throw(BioLCCCException) {
		std::vector<ChemicalGroup> parsedSequence = parseSequence(sequence, chemBasis);

		double peptideMonoisotopicMass = 0;

		for (std::vector<ChemicalGroup>::const_iterator i = parsedSequence.begin(); i < parsedSequence.end(); i++) {
			peptideMonoisotopicMass += i->monoisotopicMass();
		}

		return peptideMonoisotopicMass;
	}

}
