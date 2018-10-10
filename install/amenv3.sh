#!/usr/bin/env bash

conda env create -f `dirname "$0"`/envs/amenv3.yml
source activate amenv3
cd $(mktemp -d)
git clone --recursive https://github.com/PACChem/x2z .
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
make install
