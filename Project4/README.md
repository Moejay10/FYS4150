# Project 4 - Statistical Physics & The Ising Model

## Structure
- Codes
  - All of the Python and C++ scripts are contained in this folder
- Plots
  - All results are contained in this folder
- Tables
  - All text files are contained in this folder


## Packages
The required packages/libraries/modules to run the C++ scripts are
- Armadillo
- LAPACK
- LBLAS

The required packages to run the Python scripts are
- numpy
- matplotlib
- os


## How to run the files
Compile `main.cpp` as `c++ -O3 -fopenmp -o <executable> main.cpp Functions.cpp` and run with `./<executable>`

All plots are then produced by running `python3 plotter.py`

Unit test: Compile `unit_test.cpp` as `c++ -o <executable> test.cpp Functions.cpp` and run with `./<executable>`.
