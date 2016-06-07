#architecture file for Andromede cluster

# Debugging options
#DEBUG=-g
#DEBUG= -g -O0 -traceback -check all -check bounds -check uninit -ftrapuv - debug all -gen-interfaces -warn interfaces
DEBUG=

#flags for openmp
OPENMP = -openmp

#Choose compiler
FC=mpiifort

#ifort
F90FLAGS= -O3 -real-size 64 $(OPENMP) $(DEBUG) 

