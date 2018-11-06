#Makefile for the JAC cluster QuaLiKiz TCI interface

## Set the environment from the top-level Makefile of JETTO
## --------------------------------------------------------
include ../../include.mk

## Set the local environment for ITCequ
## ------------------------------------
include ../include.mk

# Fortran files for QUALIKIZ
OBJS_QUALIKIZ = kind.o mod_io_management.o datmat.o datcal.o dispfuncs.o asymmetry.o mod_make_io.o \
FLRterms.o mod_fluidsol.o mod_contour.o trapints.o passints.o calltrapints.o callpassints.o \
mod_fonct.o callpassQLints.o calltrapQLints.o QLflux.o nanfilter.o calcroutines.o mod_saturation.o qualikiz.o

# Math routine files
OBJS_ROUT = hrmsym.o dqagse_qlk.o cpqr79.o rpqr79.o dsort.o \
ceie.o ceik.o davint.o calci0.o calci1.o slatecroutines.o 

# CUBART routine files
OBJS_CUBART =  buckley.o internal_types.o ds_routines.o divide.o \
rule_tn.o rule_t3.o rule_t2.o rule_c2.o rule_c3.o rule_cn.o \
rule_general.o rule_1.o region_processor.o volume.o \
check.o global_all.o error_handling.o cui.o                       

FOBJ=$(OBJS_ROUT) $(OBJS_CUBART) $(OBJS_QUALIKIZ) qlk_tci_module.o

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
	ar vr $(LIBNAME) $?

# make 'realclean' option
# -----------------------
realclean:
	rm -f *.o *.mod $(LIBNAME)

# make 'clean' option
# -------------------
clean:
	rm -f *.o *.mod

# dependencies
rule_1.o: buckley.o
rule_c2.o: buckley.o
rule_t2.o: buckley.o
rule_c3.o: buckley.o
rule_t3.o: buckley.o
rule_tn.o: buckley.o
rule_cn.o: buckley.o
volume.o: buckley.o
divide.o: buckley.o internal_types.o
rule_general.o: rule_t3.o rule_t2.o rule_c2.o rule_c3.o rule_1.o \
                internal_types.o
region_processor.o : divide.o rule_general.o
global_all.o: region_processor.o volume.o ds_routines.o internal_types.o
cui.o: global_all.o check.o error_handling.o ds_routines.o internal_types.o
internal_types.o: buckley.o

%.o %.mod: %.f90
	$(FC) $(OPENMP) $(F90FLAGS) -c $<

%.o %.mod: %.f
	$(FC) $(OPENMP) $(F90FLAGS) -c $<
