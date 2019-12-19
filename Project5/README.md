# Project 5 - Partial Differential Equations

## Structure
- Codes
  - All of the Python and C++ scripts are contained in this folder
- Plots
  - All results are contained in this folder

## Packages
The required packages/libraries/modules to run the C++ scripts are
- Armadillo
- LAPACK
- LBLAS

The required packages to run the Python scripts are
- numpy
- matplotlib
- mpl_toolkits
- tqdm
- os

Compile `main.cpp` as `c++ -O3 -fopenmp -o <executable> main.cpp Functions.cpp` and run with `./<executable>`

All plots are then produced by running `python3 plotter.py`

Unit test: Compile `Unit_test.cpp` as `c++ -fopenmp -O3 -o <executable> Unit_test.cpp Functions.cpp` and run with `./<executable>`.

The folders and structure in this repository is needed for the programs to run correctly.
