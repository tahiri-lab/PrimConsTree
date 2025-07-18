#!/bin/bash

# Get FACT2 source code
git clone https://github.com/Mesh89/FACT2.git
cd FACT2/src/

# Get BOOST 1.58
wget https://archives.boost.io/release/1.58.0/source/boost_1_58_0.tar.gz
tar -xzf boost_1_58_0.tar.gz

# Compile
g++ -std=c++11 -include cstdint -I ./boost_1_58_0/ -O3 -o FACT++ Tree.cpp main.cpp
mv FACT++ ../../fact2
cd ../../

# Remove source code
rm -rf FACT2


#####
# FACT
#####

git clone https://github.com/Mesh89/FACT.git
cd FACT/src/
make
mv a ../../fact
cd ../../
rm -rf FACT
