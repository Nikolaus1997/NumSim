#!/bin/bash
sudo rm -r build/*
cd build/
cmake ..
make -j
sudo make install

echo "======================================================================="
echo "scenario: $1"
time mpirun -n 4  ./src/numsim ../$1