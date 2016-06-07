#!/bin/csh
module unload PrgEnv-cray 
module load PrgEnv-intel
cp mkfiles/Makefile_edison Makefile
cp mathlib/mkfiles/Makefile_edison mathlib/Makefile
cp cubpack/mkfiles/Makefile_edison cubpack/Makefile
make
