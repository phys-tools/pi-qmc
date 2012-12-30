#!/bin/bash

git clone git@github.com:phys-tools/pi-qmc.git

mkdir pibuild
cd pibuild
export CC=gcc-mp-4.7
export CXX=g++-mp-4.7
export CXXFLAGS="-O3 -g -Wall -ffast-math -march=native"
export CXXFLAGS+=" -ftree-vectorize -fomit-frame-pointer -pipe"
cmake ../pi-qmc
make -j2

echo; echo "  Running unit tests..."; echo
make -j2 unittest
bin/unittest_pi

cd test/system
echo; echo "Running system integration tests..."; echo
nosetests-2.7 -v --rednose --with-xunit
cd ../..


