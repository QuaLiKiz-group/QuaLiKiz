#!/bin/make -f
QUALIKIZ?=..
SRC?=$(QUALIKIZ)/src
include $(SRC)/Makefile.inc
# Directory for math routines
LIBSRC?=$(QUALIKIZ)/lib/src
CUBPACK_DIR?=$(LIBSRC)/cubpack
GENZ_DIR?=$(LIBSRC)/genz
NAG_DIR?=$(LIBSRC)/nag
SLATEC_DIR?=$(LIBSRC)/slatec
SPECFUN_DIR?=$(LIBSRC)/specfun
FUKUSHIMA_DIR?=$(LIBSRC)/fukushima

CUBPACK_INC?=-I$(CUBPACK_DIR)
CUBPACK?=$(CUBPACK_INC) -L$(CUBPACK_DIR) -lcubpack
GENZ?=-L$(GENZ_DIR) -lgenz
NAG?=-L$(NAG_DIR) -lnag
SLATEC?=-L$(SLATEC_DIR) -lslatec
SPECFUN?=-L$(SPECFUN_DIR) -lspecfun
FUKUSHIMA?=-L$(FUKUSHIMA_DIR) -lfukushima

LIBS=$(CUBPACK) $(GENZ) $(NAG) $(SLATEC) $(SPECFUN) $(FUKUSHIMA)
# Parallel or serial compilation
ifneq ($(MPI_FLAGS),-DMPI)
      OBJS_QUALIKIZ += mpivoid.o
endif

#all: $(CUB)/libcubpack.a $(ROUT)/librout.a QuaLiKiz.exe
# qlk_makeflux.exe qlk_redoQL.exe

#qlk_redoQL.exe: qlk_redoQL.f90 $(OBJS_MAKEFLUX) $(ROUT)/librout.a
#	$(FC_PREAMBLE) $(FC_WRAPPER) -o qlk_redoQL.exe qlk_redoQL.f90 $(OBJS_QUALIKIZ) ${QFLAGS} ${OPENMP} $(MPI) -L$(ROUT) -lrout -L$(CUB) -lcubpack 

#qlk_makeflux.exe: qlk_makeflux.f90 $(OBJS_MAKEFLUX) $(ROUT)/librout.a
#	$(FC_PREAMBLE) $(FC_WRAPPER) -o qlk_makeflux.exe qlk_makeflux.f90 $(OBJS_QUALIKIZ) ${QFLAGS} ${OPENMP} $(MPI) -L$(ROUT) -lrout -L$(CUB) -lcubpack 
SRCS_QUALIKIZ=calcroutines.f90 callpassQLints.f90 calltrapQLints.f90 mod_fonct.f90 QLflux.f90 qualikiz.f90 diskio.f90
OBJS_QUALIKIZ=$(SRCS_QUALIKIZ:%f90=%o)
MODS_QUALIKIZ=$(SRCS_QUALIKIZ:%f90=%mod)

SRCS_CORE=asymmetry.f90 callpassints.f90 calltrapints.f90 datcal.f90 datmat.f90 \
		  dispfuncs.f90 FLRterms.f90 kind.f90 mod_contour.f90 mod_fluidsol.f90 \
		  mod_make_io.f90 mod_saturation.f90 nanfilter.f90 \
		  passints.f90 trapints.f90
OBJS_CORE=$(SRCS_CORE:%f90=%o)
MODS_CORE=$(SRCS_CORE:%f90=%mod)

all: QuaLiKiz

# Qualikiz objects
qlflux.mod: QLflux.mod
callpassqlints.mod: callpassQLints.mod
calltrapqlints.mod: calltrapQLints.mod
qualikiz.o: mod_make_io.mod calcroutines.mod mod_saturation.mod
callpassQLints.o: callpassints.mod $(CUBPACK_DIR)/cui.mod
calltrapQLints.o: calltrapints.mod $(CUBPACK_DIR)/cui.mod
QLflux.o: kind.mod datmat.mod datcal.mod callpassqlints.mod calltrapqlints.mod
qlk_standalone.o: kind.mod
calcroutines.o: mod_fonct.mod qlflux.mod flrterms.mod mod_fluidsol.mod mod_contour.mod mod_make_io.mod asymmetry.mod nanfilter.mod
mod_fonct.o: callpassints.mod $(CUBPACK_DIR)/cui.mod calltrapints.mod
# Core objects
flrterms.mod: FLRterms.mod
callpassints.o: passints.mod 
calltrapints.o: trapints.mod
passints.o: dispfuncs.mod
trapints.o: dispfuncs.mod
dispfuncs.o: kind.mod datcal.mod datmat.mod
FLRterms.o: kind.mod datcal.mod datmat.mod
mod_fluidsol.o: kind.mod datcal.mod datmat.mod
mod_contour.o: kind.mod datcal.mod
# Makeflux objects
asymmetry.o: kind.mod datmat.mod datcal.mod
datcal.o: kind.mod
datmat.o: kind.mod
mod_make_io.o: kind.mod datmat.mod datcal.mod
mod_saturation.o: kind.mod datmat.mod datcal.mod

###############################################################################
QuaLiKiz: qlk_standalone.f90 $(OBJS_CORE) $(OBJS_QUALIKIZ) 
	$(FC_PREAMBLE) $(FC_WRAPPER) $(FFLAGS) $(MPI_FLAGS) $(OPENMP_FLAGS) $(OBJS_CORE) $(OBJS_QUALIKIZ) $(LIBS) qlk_standalone.f90 -o $@


objs_core: $(OBJS_CORE)
objs_qualikiz: $(OBJS_QUALIKIZ) objs_core
$(OBJS_QUALIKIZ) $(OBJS_CORE):%.o:%.f90
	$(FC_PREAMBLE) $(FC_WRAPPER) $(FFLAGS) $(MPI_FLAGS) $(OPENMP_FLAGS) $(LIBS) -c $<


$(MODS_QUALIKIZ) $(MODS_CORE):%.mod:%.o


clean:
	rm -f *.o *.mod


distclean: clean
	rm -f QuaLiKiz


dump_variables:
	@echo FC_PREAMBLE=$(FC_PREAMBLE)
	@echo FC_WRAPPER=$(FC_WRAPPER)
	@echo FFLAGS=$(FFLAGS)
	@echo MPI_FLAGS=$(MPI_FLAGS)
	@echo OPENMP_FLAGS=$(OPENMP_FLAGS)
	@echo
	@echo QUALIKIZ=$(QUALIKIZ)
	@echo LIBS=$(LIBS)
	@echo
	@echo MODS_CORE=$(MODS_CORE)
	@echo MODS_QUALIKIZ=$(MODS_QUALIKIZ)