MD Simulation Trajectory Analysis Tools
=====================

This git repository contains multiple useful analysis tools
that can read and analize trajectories obtained from MD simulations.

Building The Analysis Tools
===================

This project uses [trajmod](https://github.com/scychon/trajmod.git) repository for its build system. To build it, follow these steps:

1. Follow the steps to install the required [trajmod](https://github.com/scychon/trajmod.git) repository.

2. Clone or download this repository to a directory you want to build this toolkit.

3. Change the first line of install.sh file to reflect the location of trajmod module that you installed.

    export TRAJMOD=YOUR_TRAJMOD_DIRECTORY


4. Run 'bash install.sh' to compile entire toolkit.
 Alternatively, you can export the trajmod environmental variables
 (the first three lines in install.sh) and selectively run 'compile.sh'
 or 'compile_omp.sh' file in the specific analysis tool that you want to install

5. Change the structure of 'param.dat' and 'topol.dat' files to use with your own trajectory.

6. Run the analysis with your trajectory.


Test Cases
==========

The 'data' folder contains the example trajectory file that can be used to test
 tools in 'ionicLiquid' and 'equipartition' modules.
To test each tools, simply run
 
    ./TOOLNAME -f param.dat 

This should generate xxx.dat files related to the analysis module


