#!/bin/make -f
QUALIKIZ=
FC=mpif90
FFLAGS= -Mdalign -r8 -Mpreprocess -m64 -mcmodel=medium -O3
OPENMP_FLAGS=-mp=nonuma
MPI_FLAGS=-DMPI
