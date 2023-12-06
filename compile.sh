#!/bin/bash
rm -r build/*
cd build/
cmake ..
make -j
make install

echo "======================================================================="
echo "scenario: $1"
time ./src/numsim ../$1