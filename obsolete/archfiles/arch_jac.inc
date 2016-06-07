#architecture file for JAC cluster

# Debugging options
DEBUG=
#DEBUG= -g 

#pgf	
OPENMP=-mp=nonuma

# FORTRAN-90 
FC=mpif90

#pgf
F90FLAGS= -Mdalign -r8 -Mpreprocess -m64 -mcmodel=medium -O3 $(DEBUG) $(OPENMP)

