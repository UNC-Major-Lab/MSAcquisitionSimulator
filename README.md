#**MSAcquisitionSimulator**
*MSAcquisitionSimulator* is a collection of three command-line tools that simulate data-dependent acquisition algorithms on *in silico* generated ground truth data of liquid chromatography-mass spectrometry proteomics experiments. 

##**Installation Procedures**   
  
**Prerequisites**  
1. CMake 2.8+  
2. g++ 4.7+  
3. Boost 1.48+ including the compiled program_options library  
**Build**  
1. Clone the project onto your computer  
2. Navigate to the project root directory  
3. Execute the following commands:  
```ShellSession
$ cd build/
$ cmake ../src
$ make
$ cd ../bin/
$ ls
AcquisitionSimulator  FASTASampler  GroundTruthSimulator
```
if g++ and/or boost are not in your path, use this command:  
```ShellSession
$ cmake ../src -DCMAKE_CXX_COMPILER=/path/to/g++ -DBOOST_ROOT=/path/to/boost/root/boost-1.48.0/
```
The three compiled binary files, *FASTASampler*, *GroundTruthSimulator*, and *AcquisitionSimulator* should be in [project_root]/bin/  

##**Usage**  

*MSAcquisitionSimulator* comprises of three separate programs.  

*FASTASampler* takes a .fasta file as input (it can be gzipped) and outputs another .fasta file. The output contains a random subset of proteins, with each protein header appended with a '#' followed by a value for that protein's abundance. The distribution of abundance and the size of the protein subset is determined by the user.

```ShellSession
$ ./FASTASampler --help
USAGE: FASTASampler [options] input.fasta

Options:
  --help                                Print usage and exit.
  -d [ --distribution ] arg             Choose abundance distribution. Options:
                                        normal, lognormal
  -m [ --mean ] arg (=10)               The mean of the normal distribution, or
                                        m parameter for log-normal (log10).
  -s [ --stdev ] arg (=0.90000000000000002)
                                        The standard deviation of the normal 
                                        distribution, or s parameter for 
                                        log-normal (log10). 99.97% of the data 
                                        will be within +-3 standard deviations
  -p [ --percentage ] arg (=1)          Percentage of proteins to sample from. 
                                        Takes precedence over --numprot.
  -n [ --numprot ] arg                  Number of proteins to sample from 
                                        FASTA.
  -o [ --fasta_out ] arg (=sample.sim.fasta)
                                        output path for sampled FASTA file.
```

*GroundTruthSimulator* takes as input the sampled .fasta file from *FASTASampler* and a configuration file. An example configuration file is found in the root directory: *ground_truth.conf*. This program outputs the ground truth data necessary for acquisition simulation. It models post-translational modifications (PTMs), digestion, chromatographic separation, electrospray ionization, and isotopic distributions.

```ShellSession
$ ./GroundTruthSimulator --help

USAGE: GroundTruthSimulator [options] input.fasta

Options:
  --help                                Print usage and exit.
  -c [ --config ] arg (=ground_truth.conf)
                                        Input path to config file.
  -o [ --ground_truth_out_path ] arg (=ground_truth.tab)
                                        Output path for ground truth file.

```

*AcquisitionSimulator* takes as input the ground truth file from *GroundTruthSimulator* and a configuration file. An example configuration file is found in the root directory: *acquisition.conf*. This program simulates a user selected data-dependent acquisition algorithm on the ground truth data. It models chromatographic elution shape, ion accumulation, MS1 spectra, scan time durations, and database search peptide-spectrum-matches (PSMs). It currently does not simulate MS2 fragmentation spectra. The output includes an mzML file and a .fido file. The .fido file includes all PSMs in the format required for Fido to perform protein inference. Fido is available here: http://noble.gs.washington.edu/proj/fido/

```ShellSession
$ ./AcquisitionSimulator --help

USAGE: AcquisitionSimulator [options] ground_truth.tab

Options:
  --help                                Print usage and exit.
  -c [ --conf ] arg (=acquisition.conf) Input path to config file.
  -o [ --mzml_out_path ] arg (=sample.mzML)
                                        output path for mzML file.
  -f [ --fido_out_path ] arg (=sample.fido)
                                        output path for fido file.
``` 
**Example**

This is an example using 1% of the proteome:
```ShellSession
$ ./FASTASampler -dlognormal -m10 -s0.9 -p0.01 -o sampled_human_swissprot.fasta ~/Downloads/uniprot_homo_sapiens_proteome.fasta.gz 

Number of proteins in FASTA: 91618
Number of proteins sampled: 916

Abundance distribution statistics
Min: 7.74848e+06	Max: 7.40865e+12	Dynamic range: 5.98052
Median: 1.18331e+10	Mean: 8.63132e+10	Stdev: 3.98826e+11
```

```ShellSession
$ ./GroundTruthSimulator -c ../ground_truth.conf -o ground_truth_human.tab sampled_human_swissprot.fasta 

Modification registered: Carbamidomethyl car C H3,C2,N1,O1 1 0.77 0.9999 
Modification registered: Oxidation ox M O1 0.5 1.8215 0.05 
Enzyme registered: TrypsinR R P
Enzyme registered: TrypsinK K P
Enzyme registered: Random * 

Number of proteins processed: 916 of 916. Currently processing: sp|Q9NQH7|XPP3_HUMAN. Abundance: 7.63723e+10. Sequence length: 507
Number of peptides processed: 2947632 of 2947632. Number of ions passing abundance thresholds: 3552781
Elapsed time: 640 seconds

Ion abundance distribution
		x-axis: log2(count)
14	*
13	***
12	********
11	***********
10	**************
9	***************
8	*****************
7	******************
6	*******************
5	********************
4	*********************
y-axis: log10(ion abundance)

Most abundant ion for protein
		x-axis: log2(count)
14	*
13	
12	*****
11	********
10	*********
9	*********
8	*******
7	****
6	***
5	**
y-axis: log10(ion abundance)
Sorting ions by retention time. This might take a while..
```
The next command was run with the default acquisition.conf file.
```ShellSession
$ ./AcquisitionSimulator -c ../acquisition.conf -f human.fido ground_truth_human.tab

PTM registered: Carbamidomethyl C H3,C2,N1,O1
Enzyme registered: Trypsin RK P
Parsing FASTA file and digesting proteins...
Digestion complete.
Simulating Acquisition:
Current time: 10800 seconds. MS1 count: 2985. MS2 count: 26757. Num PSMs >= 0.9: 2365
Simulation Complete.
Elapsed time: 161 seconds
```
If you've installed Fido (http://noble.gs.washington.edu/proj/fido/) then you can perform protein inference on the results: 
```ShellSession
$ Fido human.fido .01 .1 .01 > human_fido_results.txt
$ less human_fido_results.txt
1 { sp|P57075|UBS3A_HUMAN }
1 { sp|Q96N96|SPT13_HUMAN }
1 { sp|Q86UQ4|ABCAD_HUMAN }
1 { sp|P07202|PERT_HUMAN }
1 { sp|Q9Y5L0|TNPO3_HUMAN }
1 { sp|Q8IX01|SUGP2_HUMAN }
1 { sp|Q9P2G1|AKIB1_HUMAN }
1 { sp|Q96RG2|PASK_HUMAN }
...
$ less human_fido_results.txt | grep "^1" | wc -l
68
$ less human_fido_results.txt | grep "^0.9" | wc -l
113
```

##**Simulator Details**  
###**Ground truth generation**  
**Digestion**  
*In silico* digestion is not performed in the same way as it is in a typical database search engine. Instead of having hard cut-offs for the maximum number of missed cleavages and minimum number of enzymatic termini, *GroundTruthSimulator* calculates the probability of each peptide's existence based on each enzyme's probability of cleavage. For a particular peptide, the probability of its existence from a single copy of a protein = 1 - PROB(not existing) = 1 - ( PROB(no cleavage at n-term) * PROB(no cleavage at c-term) * PRODUCT_over_all_missed_cleavages(PROB(cleavage)) ). To determine the number of copies that will exist of this peptide, we multiply the final digestion probability by the protein abundance.  
**Post-translational modifications**  
Similarly to digestion, *GroundTruthSimulator* does not have hard cut-offs for the maximum number of dynamic modifications. Even the idea of dynamic vs static modifications is not used. Each modification has a user-defined probability of occupying a particular site (e.g. there are probably serines that are *never* phosphorylated), and a user-defined percentage of that residue that will be modified (e.g. if a particular serine *is* chosen to be phosphorylated, maybe only 1% of it ever exists in that form at any given time). A "static" modification can be simulated by setting both values to 100%. Candidate modification sites are first randomly assigned to each protein based on each PTM's probability of occupying a particular site. For each peptide created during the digestion process, modification combinations are created and their probability existing from a single protein is calculated. Combinations with low abundance are pruned.  The probability of a particular modification combination = PRODUCT_over_all_modification_states(PROB(modification state)). The probability of the modification state being "no modification" for a specific peptide's residue is 1 - SUM_over_all_modifications(PROB(modification)).  
**Retention time**  
Determined via the BioLCCC library (http://pythonhosted.org/pyteomics.biolccc/)  
**Elution shape**  
The shape of an ion's elution profile is modeled by an Exponential Gaussian Hybrid (EGH) function as described in “A Hybrid of Exponential and Gaussian Functions as a Model of Asymmetric Chromatographic Peaks”, by Kevin Lan and James W. Jorgenson, Journal of Chromatography, Part A, 915, 1-13 (2001).  
**Ionization efficiency**  
The probability of a peptide ionizing is sampled from a uniform random distribution between 0 and 1.  
**Charge state distribution**  
The probability a peptide is has a charge of k is equal to a binomial distribution with n = the number of basic residues + 1 (for the n-terminus) and probability of success p = .7 + .3 x uniform_random(0,1).  
**Isotopic distribution**  
An ion's isotopic distribution is determined with the libmercury++ library based on Rockwood, A.L. and Haimi, P.: "Efficent calculation of Accurate Masses of Isotopic Peaks", Journal of The American Society for Mass Spectrometry. JASMS 03-2263, 2006.  
**Ion abundance**  
We assume that digestion, modifications, ionization efficiency, charge, and isotopic distributions are all independent of each other. Therefore, an ion with charge z and isotope m has total abundance = Protein abundance x prob(digestion) x prob(PTMs) x prob(ionization efficiency) x prob(charge state = z) x prob(isotope = m)   
If at any point in the simulation, the abundance gets below the filter threshold set in ground_truth.conf, then it is filtered and not included in the latter stages of the simulator.
###**Acquisition simulation**
####**MS1 Scan**  
**Ion abundance**  
The elution shape is numerically integrated using Simpson's method one millisecond at a time for every ion present at the current time and m/z constraints. This integration continues until we've reached the target total ion count, or the maximum injection time.  
**Scan time**  
The elapsed time for a scan is equal to max(injection_time, transient_time) + scan_overhead_time. This models the scan time for a QExactive-like instrument.  
**Raw signals**  
The Cauchy-Lorentz distribution is used to model the peak shape for each ion. The raw signal is a mixture of these distributions and therefore the intensity at each m/z is equal to the sum of the contribution from each ion's distribution. These signals are then centroided.  
####**MS2 Scan**  
Raw signals are generated for the precursor ions only. Fragmentation is not modeled. Scan time and ion abundance are computed identically as an MS1 scan.  
**Sequence determination**  
First, the precursor ion fraction (PIF) is calculated for each peptide in the scan. The PIF of a peptide is defined as the sum of ion intensities for all ions of that peptide (i.e. the sum of all isotope intensities for that peptide) divided by the total ion intensity of the scan. Next, a peptide is randomly selected from these peptides - weighted by their PIF. If the peptide is in the peptide database used to simulate a database search, and within the user-defined mass tolerance, then the peptide-spectrum-match (PSM) sequence is set to this peptide. Otherwise, we randomly choose if the PSM maps to a decoy (50% probability). If it's a decoy, the peptide sequence is "DECOY_#". Otherwise, we randomly (uniformly) choose a sequence from the peptides in the database search that are within our mass tolerance of the targeted precursor m/z.  
**Probability determination**  
 If the PSM sequence was mapped to a peptide in the scan, then the PSM probability is set to that peptide's PIF. Otherwise, the PSM probability is sampled from a truncated exponential distribution between 0 and 1. This null distribution's lambda parameter is set by the user in the configuration file. Reasonable defaults for all parameters are provided.  
####**Acquisition loop**  
The acquisition simulator can be described by this pseudocode (nearly identical to the actual code). The controller (an acquisition algorithm) gives scan requests, the oracle generates data for that request, the controller processes that data, and the current time is updated.  
```cpp
while (current_time < acquisition_length) {
	scanRequest = controller.get_scan_request();
	scan = oracle.get_scan_data(scan_request, current_time);
	controller.process_scan(scan);
	current_time += scan.elapsed_time;
}
```
##**Creating a Custom Acquisition Controller**

* Create a class that inherits the *AcquisitionController* interface.  

* Implement a constructor that takes in a *std::vector\<std::string\>* parameter which should be parsed by boost::program_options. Use *AbstractTopN* as an example.  

* Implement a function to create a new scan request (std::unique_ptr\<ScanRequest\> get_scan_request(double current_time))  

* A ScanRequest consists of the min_mz and max_mz for an MS1 and also a BasicPeak (mz, intensity) for the target precursor peak, and parent_scan_id for an MS2. Again see *AbstractTopN* for an example.  

* Implement a function to process a new scan (void process_scan(Scan* scan)).  

* An MS1Scan consists of a *std::vector\<BasicPeak\>*, time_of_scan, elapsed_time, scan_id, and scan_type (MS1 or MS2).  

* An MS2Scan also contains the PIF for each peptide, total ion intensity, the targeted precursor peak, the PSM sequence and probability, and a list of proteins the peptide maps to.  
		
* Register your controller with the get_controller() function in mainAcquisitionSimulator.cpp so that it will be used when the acquisition.conf file provides the name of your controller.
