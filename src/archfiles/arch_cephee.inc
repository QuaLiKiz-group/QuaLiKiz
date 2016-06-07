#architecture file for Cephee cluster

# MPI OPTION
MPI=-DMPI
#MPI=

# Debugging options
#DEBUG=
#DEBUG= -g -O0 -traceback -check all -check bounds -check uninit -ftrapuv - debug all -gen-interfaces -warn interfaces

#gfortran
#DEBUG = -fbacktrace -ffpe-trap=invalid -ffpe-trap=zero -ffpe-trap=overflow -ffpe-trap=underflow -fdump-core
DEBUG=


ifeq ($(MPI),-DMPI)
 
   #mpi include directory
   #INCLUDE_MPI=-I/usr/mpi/gcc/openmpi-1.8.4/include -I/usr/mpi/gcc/openmpi-1.8.4/lib64

   #mpi library
   #pgf
   #LIB_MPI=-L/usr/mpi/gcc/openmpi-1.8.4/lib64 -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -ldl

   ##gfortran
   #LIB_MPI=-L/usr/mpi/gcc/openmpi-1.8.4/lib64 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -ldl

   #flags for openmp

   #gfortran
   OPENMP=-fopenmp 

   #pgf	
   #OPENMP=-mp=nonuma
   
   #ifort
   #OPENMP = -openmp
else
  
   INCLUDE_MPI=
   LIB_MPI=
   OPENMP=
endif


# FORTRAN-90 
#FC=pgf95
#FC=/usr/mpi/gcc/openmpi-1.8.4/bin/mpifort
#FC=ifort
FC=mpifort

#pgf
#F90FLAGS= -Mdalign -r8 -Mpreprocess -O3 $(MPI) $(INCLUDE_MPI) $(INCLUDE_QLK) $(DEF) $(DEBUG)

#gfortran
#F90FLAGS= -O3 -ffree-line-length-none -fdefault-real-8 -finit-local-zero $(MPI) $(INCLUDE_MPI) $(DEBUG)
F90FLAGS= -O3 -ffree-line-length-none -fdefault-real-8 $(MPI) $(INCLUDE_MPI) $(DEBUG)

#ifort
#F90FLAGS= -fast -no-ipo -real-size 64 $(MPI) $(INCLUDE_MPI) $(DEBUG)

F90LIBS=$(LIB_MPI)
QFLAGS = $(F90FLAGS) $(F90LIBS) -Icubpack 
