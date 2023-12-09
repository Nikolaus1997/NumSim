#!/bin/bash
sudo rm -r build/*
cd build/
cmake ..
make -j
sudo make install

echo "======================================================================="
echo "scenario: $1"
time mpirun -n 2  ./src/numsim_parallel ../$1