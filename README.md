#**MSAcquisitionSimulator**

##**Installation Procedures**   
  
**Prerequisites**  
1. CMake 2.8+  
2. g++ 4.7+  
3. Boost 1.48+ including the compiled program_options library  
**Build**  
1. Clone the project onto your computer  
2. Navigate to the project root directory  
3. Execute the following commands:  
```Shell
$ cd build/
$ cmake ../src
$ make
$ cd ../bin/
$ ls
AcquisitionSimulator  FASTASampler  GroundTruthSimulator
```
if g++ and/or boost are not in your path, use this command:  
```Shell
$ cmake ../src -DCMAKE_CXX_COMPILER=/path/to/g++ -DBOOST_ROOT=/path/to/boost/root/boost-1.48.0/
```
The three compiled binary files, *FASTASampler*, *GroundTruthSimulator*, and *AcquisitionSimulator* should be in [project_root]/bin/  
*** 
##**Usage**  

MSAcquisitionSimulator comprises of three separate programs. First, *FASTASampler* takes a .fasta file as input and outputs another .fasta file. The output contains a subset of proteins, with each protein header appended with a '#' followed by a value for that protein's abundance. The distribution of abundance and the size of the protein subset is determined by the user.

```Shell
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

##**Simulator Details**

##**Creating a Custom Acquisition Controller**
