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
***
##**Simulator Details**
***
##**Creating a Custom Acquisition Controller**
