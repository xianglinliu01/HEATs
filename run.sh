#!/bin/bash

# To use different parameters, such as unit cell size, edit ./Python/genNeighbors.py
cd Python
python genNeighbors.py
cd ..
cd C++
#mpicxx -O3 -std=c++11 MonteCarlo.cpp
make
mpirun -n 4 ./MonteCarlo > ../Output/output.dat
cd ..
