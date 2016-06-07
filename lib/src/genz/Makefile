#!/bin/make -f
QUALIKIZ?=../../..
SRC?=$(QUALIKIZ)/src
include $(SRC)/Makefile.inc

# Fortran routine files
SRCS_GENZ=hrmsym.f
OBJS_GENZ=$(SRCS_GENZ:%f=%o)

##############################################################################
libgenz.a: $(OBJS_GENZ)
	$(AR) -cr $@ $^

all: $(OBJS_GENZ)


$(OBJS_GENZ):%.o:%.f
	$(FC_PREAMBLE) $(FC_WRAPPER) $(FFLAGS) $(MPI_FLAGS) $(OPENMP_FLAGS) -c $<


clean:
	rm -f *.o


distclean: clean
	rm -f *.a
