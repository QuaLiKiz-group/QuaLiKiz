#Makefile for the JAC cluster QuaLiKiz TCI interface

## Set the environment from the top-level Makefile of JETTO
## --------------------------------------------------------
include ../../../include.mk

## Set the local environment for ITCequ
## ------------------------------------
include ../../include.mk

# Fortran files for QUALIKIZ
OBJS_QUALIKIZ = kind.o mod_io_management.o datmat.o datcal.o dispfuncs.o asymmetry.o mod_make_io.o \
FLRterms.o mod_fluidsol.o mod_contour.o trapints.o passints.o calltrapints.o callpassints.o \
mod_fonct.o callpassQLints.o calltrapQLints.o QLflux.o nanfilter.o calcroutines.o mod_saturation.o qualikiz.o

# Math routine files
OBJS_ROUT = hrmsym.o dqagse_qlk.o cpqr79.o rpqr79.o dsort.o \
ceie.o ceik.o davint.o calci0.o calci1.o slatecroutines.o 

OBJS_ROUT2 = $(addprefix mathlib/,$(OBJS_ROUT))

# CUBART routine files
OBJS_CUBART =  buckley.o internal_types.o ds_routines.o divide.o \
rule_tn.o rule_t3.o rule_t2.o rule_c2.o rule_c3.o rule_cn.o \
rule_general.o rule_1.o region_processor.o volume.o \
check.o global_all.o error_handling.o cui.o                       

OBJS_CUBART2 = $(addprefix cubpack/,$(OBJS_CUBART))

FOBJ=$(OBJS_ROUT2) $(OBJS_CUBART2) $(OBJS_QUALIKIZ) qlk_tci_module.o

## set PHONY targets
## -----------------
.PHONY: all realclean clean $(MASTER_RULES) $(LIBNAME)

## make 'all' option
## -----------------
#all:
#	($(MAKE) -C ../../ all) || exit $$?

## Set rules passed to top-level Makefile
## --------------------------------------
#$(MASTER_RULES):
#	@($(MAKE) -C ../../ $@) || exit $$?

# Rule for libtransport.a
# -----------------------
$(LIBNAME): $(FOBJ)
	ar vr $(LIBNAME) $(notdir $?)

# make 'realclean' option
# -----------------------
realclean:
	rm -f *.o *.mod $(LIBNAME)

# make 'clean' option
# -------------------
clean:
	rm -f *.o *.mod

# dependencies
cubpack/rule_1.o: cubpack/buckley.o
cubpack/rule_c2.o: cubpack/buckley.o
cubpack/rule_t2.o: cubpack/buckley.o
cubpack/rule_c3.o: cubpack/buckley.o
cubpack/rule_t3.o: cubpack/buckley.o
cubpack/rule_tn.o: cubpack/buckley.o
cubpack/rule_cn.o: cubpack/buckley.o
cubpack/volume.o: cubpack/buckley.o
cubpack/divide.o: cubpack/buckley.o cubpack/internal_types.o
cubpack/rule_general.o: cubpack/rule_1.o cubpack/rule_t3.o cubpack/rule_t2.o \
                        cubpack/rule_c2.o cubpack/rule_c3.o cubpack/internal_types.o
cubpack/region_processor.o : cubpack/divide.o cubpack/rule_general.o
cubpack/global_all.o: cubpack/region_processor.o cubpack/volume.o cubpack/ds_routines.o cubpack/internal_types.o
cubpack/cui.o: cubpack/global_all.o cubpack/check.o cubpack/error_handling.o cubpack/ds_routines.o cubpack/internal_types.o
cubpack/internal_types.o: cubpack/buckley.o
mod_fonct.o: cubpack/cui.o

%.o %.mod: %.f90
	$(FC) $(OPENMP) $(F90FLAGS) -c $<

%.o %.mod: %.f
	$(FC) $(OPENMP) $(F90FLAGS) -c $<
