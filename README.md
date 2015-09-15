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
**Peptide generation**  
**Post-translational modifications**  
**Digestion**  
**Retention time**  
**Ionization efficiency**  
**Charge state distribution**  
**Isotopic distribution**  
**Ion abundance**  
###**Acquisition simulation**
**MS1 Scan**  
**Ion abundance**  
**Raw signals**  
**Scan time**  
**MS2 Scan**  
**Sequence determination**  
**Probability determination**  
**Acquisition loop**  

##**Creating a Custom Acquisition Controller**
