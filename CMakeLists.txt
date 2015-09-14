cmake_minimum_required(VERSION 3.2)
project(MSAcquisitionSimulator)

# The version number.
set(MSAcquisitionSimulator_VERSION_MAJOR 1)
set(MSAcquisitionSimulator_VERSION_MINOR 0)

set(FASTASampler_VERSION_MAJOR 1)
set(FASTASampler_VERSION_MINOR 0)

set(GroundTruthSimulator_VERSION_MAJOR 1)
set(GroundTruthSimulator_VERSION_MINOR 0)

set(AcquisitionSimulator_VERSION_MAJOR 1)
set(AcquisitionSimulator_VERSION_MINOR 0)

# Configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${PROJECT_SOURCE_DIR}/MSAcquisitionSimulatorConfig.h.in"
        "${PROJECT_BINARY_DIR}/MSAcquisitionSimulatorConfig.h"
)
include_directories("${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -lboost_program_options")

set(SOURCE_FILES_BIOLCCC BioLCCC/biolccc.h BioLCCC/auxiliary.h BioLCCC/chemicalbasis.h BioLCCC/chromoconditions.h BioLCCC/gradient.h BioLCCC/chemicalgroup.h BioLCCC/chain_model.h BioLCCC/gradientpoint.h BioLCCC/parsing.h BioLCCC/rod_model.h BioLCCC/biolcccexception.h BioLCCC/chemicalbasis.cpp BioLCCC/biolccc.cpp BioLCCC/gradient.cpp BioLCCC/chromoconditions.cpp BioLCCC/chemicalgroup.cpp BioLCCC/chain_model.cpp BioLCCC/gradientpoint.cpp BioLCCC/rod_model.cpp BioLCCC/parsing.cpp BioLCCC/biolcccexception.cpp BioLCCC/auxiliary.cpp)
set(SOURCE_FILES_SAMPLER mainFASTASampler.cpp Protein.cpp Protein.h FASTAParser.cpp FASTAParser.h PTM.cpp PTM.h MolecularFormula.cpp MolecularFormula.h)
set(SOURCE_FILES_GROUND_TRUTH mainGroundTruthSimulator.cpp GroundTruthText.cpp GroundTruthText.h IonDAO.cpp IonDAO.h Protein.cpp Protein.h FASTAParser.cpp FASTAParser.h PTM.cpp PTM.h MolecularFormula.cpp MolecularFormula.h Globals.cpp Globals.h Element.cpp Element.h Residue.cpp Residue.h SubatomicParticle.cpp SubatomicParticle.h Enzyme.cpp Enzyme.h libmercury++.h libmercury++.cpp Molecule.cpp Molecule.h DefaultPTM.cpp DefaultPTM.h DefaultEnzyme.cpp DefaultEnzyme.h RandomFactory.cpp RandomFactory.h PeptideGenerator.cpp PeptideGenerator.h ModificationSimulator.cpp ModificationSimulator.h DigestionSimulator.cpp DigestionSimulator.h Peptide.cpp Peptide.h PeptideForDigestion.cpp PeptideForDigestion.h JointMultinomialEnumerator.cpp JointMultinomialEnumerator.h PTMLocation.cpp PTMLocation.h PeptideDigested.cpp PeptideDigested.h IonizationSimulator.cpp IonizationSimulator.h Histogram.cpp Histogram.h ElutionShapeSimulator.cpp ElutionShapeSimulator.h)
set(SOURCE_FILES_ACQUISITION mainAcquisitionSimulator.cpp Peptide.cpp Peptide.h Protein.cpp Protein.h FASTAParser.cpp FASTAParser.h MolecularFormula.cpp MolecularFormula.h Globals.cpp FASTAParser.cpp FASTAParser.h PTM.cpp PTM.h MS2ScanRequest.cpp MS2ScanRequest.h SubatomicParticle.cpp SubatomicParticle.h libmercury++.h libmercury++.cpp Molecule.cpp Molecule.h TopNParameters.cpp TopNParameters.h StochasticN.cpp StochasticN.h ExclusionList.cpp ExclusionList.h Enzyme.cpp Enzyme.h DynamicExclusionList.cpp DynamicExclusionList.h StaticExclusionList.cpp StaticExclusionList.h Globals.h Instrument.cpp Instrument.h AcquisitionController.cpp AcquisitionController.h Scan.cpp Scan.h ScanRequest.cpp ScanRequest.h Oracle.cpp Oracle.h Sequencer.cpp Sequencer.h ElutionShapeSimulator.cpp ElutionShapeSimulator.h TopN.cpp TopN.h MS1Scan.cpp MS1Scan.h MS2Scan.cpp MS2Scan.h Peak.cpp Peak.h SpectrumAnalyzer.cpp SpectrumAnalyzer.h StaticExclusionListEntry.cpp StaticExclusionListEntry.h DynamicExclusionListEntry.cpp DynamicExclusionListEntry.h StaticExclusionMatchQuery.cpp StaticExclusionMatchQuery.h DynamicExclusionMatchQuery.cpp DynamicExclusionMatchQuery.h AbstractTopN.cpp AbstractTopN.h RandomN.cpp RandomN.h WeightedShuffle.cpp WeightedShuffle.h BasicPeak.cpp BasicPeak.h MzMLWriter.cpp MzMLWriter.h base64.cpp base64.h DBEnzyme.cpp DBEnzyme.h Enumerator.cpp Enumerator.h Element.h Element.cpp Residue.h Residue.cpp DBPTM.cpp DBPTM.h MS2Peptide.cpp MS2Peptide.h FidoWriter.cpp FidoWriter.h GroundTruthText.cpp GroundTruthText.h IonDAO.cpp IonDAO.h)

add_executable(FASTASampler ${SOURCE_FILES_SAMPLER})
add_executable(GroundTruthSimulator ${SOURCE_FILES_GROUND_TRUTH} ${SOURCE_FILES_BIOLCCC})
add_executable(AcquisitionSimulator ${SOURCE_FILES_ACQUISITION} )

find_package(ZLIB)
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    target_link_libraries(FASTASampler ${ZLIB_LIBRARIES})
    target_link_libraries(GroundTruthSimulator ${ZLIB_LIBRARIES})
    target_link_libraries(AcquisitionSimulator ${ZLIB_LIBRARIES})
endif()

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)

#find_package(Boost)
find_package(Boost COMPONENTS program_options)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    target_link_libraries(FASTASampler ${Boost_LIBRARIES})
    target_link_libraries(GroundTruthSimulator ${Boost_LIBRARIES})
    target_link_libraries(AcquisitionSimulator ${Boost_LIBRARIES})
endif()

# add the install targets
install (TARGETS FASTASampler GroundTruthSimulator AcquisitionSimulator DESTINATION bin)


