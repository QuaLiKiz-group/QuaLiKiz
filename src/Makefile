#!/bin/make -f
QUALIKIZ?=..
QUALIKIZ_SRC?=$(QUALIKIZ)/src
include $(QUALIKIZ_SRC)/Makefile.inc
# Directory for math routines
QUALIKIZ_LIBSRC?=$(QUALIKIZ)/lib/src
SLATEC_DIR?=$(QUALIKIZ_LIBSRC)/slatec
SPECFUN_DIR?=$(QUALIKIZ_LIBSRC)/specfun
FUKUSHIMA_DIR?=$(QUALIKIZ_LIBSRC)/fukushima

SLATEC?=-L$(SLATEC_DIR) -lslatec
SPECFUN?=-L$(SPECFUN_DIR) -lspecfun
FUKUSHIMA?=-L$(FUKUSHIMA_DIR) -lfukushima

QUALIKIZ_LIBS=$(SLATEC) $(SPECFUN) $(FUKUSHIMA) $(CUBA)
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

QUALIKIZ_SRCS_CORE=asymmetry.f90 callpassints.f90 calltrapints.f90 calltotalints.f90 datcal.f90 datmat.f90 \
		  dispfuncs.f90 hcub.f90 clencurt.f90 pcub.f90 FLRterms.f90 kind.f90 mod_contour.f90 mod_fluidsol.f90 \
		  mod_make_io.f90 mod_saturation.f90 nanfilter.f90 \
		  passints.f90 trapints.f90
QUALIKIZ_OBJS_CORE=$(QUALIKIZ_SRCS_CORE:%f90=%o)
QUALIKIZ_MODS_CORE=$(QUALIKIZ_SRCS_CORE:%f90=%mod)

INT_DRIV_SRCS_CORE=kind.f90 diskio.f90 datcal.f90 datmat.f90 FLRterms.f90 mod_fluidsol.f90 mod_contour.f90 dispfuncs.f90 hcub.f90 clencurt.f90 pcub.f90
INT_DRIV_OBJS_CORE=$(INT_DRIV_SRCS_CORE:%f90=%o) 
INT_DRIV_MODS_CORE=$(INT_DRIV_SRCS_CORE:%f90=%mod)

# Qualikiz objects
qlflux.mod: QLflux.mod
callpassqlints.mod: callpassQLints.mod
calltrapqlints.mod: calltrapQLints.mod
qualikiz.o: mod_make_io.mod calcroutines.mod mod_saturation.mod
callpassQLints.o: callpassints.mod hcub.mod pcub.mod
calltrapQLints.o: calltrapints.mod hcub.mod pcub.mod
QLflux.o: kind.mod datmat.mod datcal.mod callpassqlints.mod calltrapqlints.mod
qlk_standalone.o: kind.mod
calcroutines.o: mod_fonct.mod qlflux.mod flrterms.mod mod_fluidsol.mod mod_contour.mod mod_make_io.mod asymmetry.mod nanfilter.mod
mod_fonct.o: callpassints.mod calltrapints.mod calltotalints.mod hcub.mod pcub.mod
# Core objects
flrterms.mod: FLRterms.mod
callpassints.o: passints.mod 
calltrapints.o: trapints.mod
calltotalints.o: passints.mod trapints.mod
passints.o: dispfuncs.mod
trapints.o: dispfuncs.mod
dispfuncs.o: kind.mod datcal.mod datmat.mod
hcub.o: kind.mod
clencurt.o : kind.mod
pcub.o: kind.mod clencurt.mod
FLRterms.o: kind.mod datcal.mod datmat.mod
mod_fluidsol.o: kind.mod datcal.mod datmat.mod hcub.mod pcub.mod
mod_contour.o: kind.mod datcal.mod
# Makeflux objects
asymmetry.o: kind.mod datmat.mod datcal.mod hcub.mod pcub.mod
datcal.o: kind.mod
datmat.o: kind.mod
mod_make_io.o: kind.mod datmat.mod datcal.mod
mod_saturation.o: kind.mod datmat.mod datcal.mod

###############################################################################
QuaLiKiz: qlk_standalone.f90 $(QUALIKIZ_OBJS_CORE) $(OBJS_QUALIKIZ) 
	$(FC_PREAMBLE) $(FC_WRAPPER) $(FFLAGS) $(MPI_FLAGS) $(OPENMP_FLAGS) $(QUALIKIZ_OBJS_CORE) $(OBJS_QUALIKIZ) $(QUALIKIZ_LIBS) qlk_standalone.f90 -o $@

objs_core: $(QUALIKIZ_OBJS_CORE)
objs_qualikiz: $(OBJS_QUALIKIZ) objs_core
$(OBJS_QUALIKIZ) $(QUALIKIZ_OBJS_CORE):%.o:%.f90
	$(FC_PREAMBLE) $(FC_WRAPPER) $(FFLAGS) $(MPI_FLAGS) $(OPENMP_FLAGS) $(QUALIKIZ_LIBS) -c $<


$(MODS_QUALIKIZ) $(QUALIKIZ_MODS_CORE):%.mod:%.o


clean:
	rm -f *.o *.mod *__genmod.f90


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
	@echo QUALIKIZ_LIBS=$(QUALIKIZ_LIBS)
	@echo
	@echo QUALIKIZ_MODS_CORE=$(QUALIKIZ_MODS_CORE)
	@echo MODS_QUALIKIZ=$(MODS_QUALIKIZ)
