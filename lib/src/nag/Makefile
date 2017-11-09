#!/bin/make -f
QUALIKIZ?=../../..
QUALIKIZ_SRC?=$(QUALIKIZ)/src
include $(QUALIKIZ_SRC)/Makefile.inc


# Fortran routine files
SRCS_NAG=d01ahf.f d01ahu.f d01ahv.f d01ahw.f d01ahx.f d01ahy.f d01ahz.f d01fcf.f \
		  g05caf.f g05cay.f g05caz.f p01abf.f p01abz.f x02ajf.f x02amf.f x02bbf.f \
		  x04aaf.f x04baf.f
OBJS_NAG=$(SRCS_NAG:%f=%o)

##############################################################################
libnag.a: $(OBJS_NAG)
	$(AR) -cr $@ $^

all: $(OBJS_NAG)


$(OBJS_NAG):%.o:%.f
	$(FC_PREAMBLE) $(FC_WRAPPER) $(FFLAGS) $(MPI_FLAGS) $(OPENMP_FLAGS) -c $<


clean:
	rm -f *.o


distclean: clean
	rm -f *.a
