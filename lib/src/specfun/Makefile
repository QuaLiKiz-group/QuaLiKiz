#!/bin/make -f
QLKDIR?=../../..
QUALIKIZ_SRC?=$(QLKDIR)/src
include $(QUALIKIZ_SRC)/Makefile.inc

# Fortran routine files
SRCS_SPECFUN=calci0.f90 calci1.f90
OBJS_SPECFUN=$(SRCS_SPECFUN:%f90=%o)

##############################################################################
libspecfun.a: $(OBJS_SPECFUN)
	$(AR) -cr $@ $^


all: $(OBJS_SPECFUN)


$(OBJS_SPECFUN):%.o:%.f90
	$(FC_PREAMBLE) $(FC_WRAPPER) $(FFLAGS) $(MPI_FLAGS) $(OPENMP_FLAGS) -c $<


clean:
	rm -f *.o *.mod *__genmod.f90


distclean: clean
	rm -f *.a
