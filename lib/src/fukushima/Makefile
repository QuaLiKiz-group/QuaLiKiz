#!/bin/make -f
QUALIKIZ?=../../..
SRC?=$(QUALIKIZ)/src
include $(SRC)/Makefile.inc

# Fortran routine files
SRCS_FUKUSHIMA_MAY=elbd.f90 relbd.f90 celbd.f90 rcelbd.f90
OBJS_FUKUSHIMA_MAY=$(SRCS_FUKUSHIMA_MAY:%f90=%o)
SRCS_FUKUSHIMA_FEB=ceik.f90 ceie.f90 ceib.f90 ceid.f90 ceis.f90
OBJS_FUKUSHIMA_FEB=$(SRCS_FUKUSHIMA_FEB:%f90=%o)

##############################################################################


libfukushima.a: $(OBJS_FUKUSHIMA_MAY) $(OBJS_FUKUSHIMA_FEB)
	$(AR) -cr $@ $^


xelbd: xelbd.f90 libfukushima.a
	$(FC_PREAMBLE) $(FC_WRAPPER) $(FFLAGS) $(MPI_FLAGS) $(OPENMP_FLAGS) xelbd.f90 libfukushima.a -o $@


xcei: xcei.f90 libfukushima.a
	$(FC_PREAMBLE) $(FC_WRAPPER) $(FFLAGS) $(MPI_FLAGS) $(OPENMP_FLAGS) xcei.f90 libfukushima.a -o $@


all: libfukushima.a xelbd xcei


$(OBJS_FUKUSHIMA_MAY) $(OBJS_FUKUSHIMA_FEB):%.o:%.f90
	$(FC_PREAMBLE) $(FC_WRAPPER) $(FFLAGS) $(MPI_FLAGS) $(OPENMP_FLAGS) -c $<


clean:
	/bin/rm -f *.o *.mod


distclean: clean
	/bin/rm -f *.a xelbd xcei
