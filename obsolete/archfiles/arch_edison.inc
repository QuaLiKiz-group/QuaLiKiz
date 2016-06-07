#architecture file for Andromede cluster

# MPI OPTION
MPI=-DMPI
#MPI=

# Debugging options
#DEBUG=-g
#DEBUG= -g -O0 -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -warn interfaces
#DEBUG = -debug all
DEBUG=

ifeq ($(MPI),-DMPI)
 
   #mpi include directory

   #INCLUDE_MPI = -I/$INTEL_HOME/mpi/intel64/lib -I/$INTEL_HOME/compiler/include 
   #LIB_MPI = -L/$INTEL_HOME/mpi/intel64/lib -lmpifort

   #INCLUDE_MPI = -I/Applications/intel2016/compilers_and_libraries_2016.1.150/linux/mpi/intel64/include -L/Applications/intel2016/compilers_and_libraries_2016.1.150/linux/mpi/intel64/lib
   INCLUDE_MPI = 

   #mpi library
   #pgf
   #LIB_MPI=-L/usr/mpi/gcc/openmpi-1.8.4/lib64 -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -ldl

   ##gfortran
   #INCLUDE_MPI=-I/usr/lib64/openmpi/bin -I/usr/lib64/openmpi/lib -I/usr/include/openmpi-x86_64
   #LIB_MPI=-L/usr/lib64/openmpi/lib -L/usr/include/openmpi-x86_64 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -ldl

   #flags for openmp

   #gfortran
   #OPENMP=-fopenmp 

   #pgf	
   #OPENMP=-mp=nonuma
   
   #ifort
   OPENMP = -openmp
else
  
   INCLUDE_MPI=
   LIB_MPI=
   OPENMP=
endif


# FORTRAN-90 
#FC=pgf90
#FC=mpif90
#FC=ifort
#FC=mpifort
#FC=mpiifort
FC=ftn

#pgf
#F90FLAGS= -Mdalign -r8 -Mpreprocess -O3 $(MPI) $(INCLUDE_MPI) $(INCLUDE_QLK) $(DEF) $(DEBUG)

#gfortran
#F90FLAGS= -O3 -ffree-line-length-none -fdefault-real-8 -finit-local-zero $(MPI) $(INCLUDE_MPI) $(DEBUG)

#ifort
F90FLAGS= -O2 -real-size 64 $(MPI) $(INCLUDE_MPI) $(DEBUG)

F90LIBS=$(LIB_MPI)
QFLAGS = $(F90FLAGS) $(F90LIBS) -Icubpack
