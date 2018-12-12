#!/bin/make -f
QLKDIR?=../../..
QUALIKIZ_SRC?=$(QLKDIR)/src
include $(QUALIKIZ_SRC)/Makefile.inc

FFLAGS_SLATEC=-fdefault-real-8  

# Fortran routine files
SRCS_SLATEC=cpqr79.f davint.f dqagse_qlk.f dsort.f rd.f rf.f rpqr79.f slatecroutines.f
OBJS_SLATEC=$(SRCS_SLATEC:%f=%o)

##############################################################################
libslatec.a: $(OBJS_SLATEC)
	$(AR) -cr $@ $^


all: $(OBJS_SLATEC)


$(OBJS_SLATEC):%.o:%.f
	$(FC_PREAMBLE) $(FC_WRAPPER) $(FFLAGS) $(MPI_FLAGS) $(OPENMP_FLAGS) -c $<


clean:
	rm -f *.o *.mod *__genmod.f90


distclean: clean
	rm -f *.a
