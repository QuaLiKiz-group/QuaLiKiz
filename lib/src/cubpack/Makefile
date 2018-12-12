#!/bin/make -f
QLKDIR=../../..
QUALIKIZ_SRC?=$(QLKDIR)/src
include $(QUALIKIZ_SRC)/Makefile.inc

SRCS_CUBPACK=buckley.f90 check.f90 cui.f90 divide.f90 ds_routines.f90 \
			 error_handling.f90 global_all.f90 internal_types.f90 region_processor.f90 \
			 rule_1.f90 rule_c2.f90 rule_c3.f90 rule_cn.f90 rule_general.f90 \
			 rule_t2.f90 rule_t3.f90 rule_tn.f90 volume.f90

OBJS_CUBPACK=$(SRCS_CUBPACK:%f90=%o)
MODS_CUBPACK=$(SRCS_CUBPACK:%f90=%mod)

##############################################################################
libcubpack.a:
# Special mod rules
extended.mod integers.mod low.mod normal.mod precision_model.mod: buckley.o
volume_computation.mod: volume.o
check_input.mod: check.o
subdivisions.mod: divide.o
quadraturerule.mod: rule_1.o
cubaturerule_general.mod: rule_general.o
cubaturerule_c2.mod: rule_c2.o
cubaturerule_c3.mod: rule_c3.o
cubaturerule_cn.mod: rule_cn.o
cubaturerule_t2.mod: rule_t2.o
cubaturerule_t3.mod: rule_t3.o
cubaturerule_tn.mod: rule_tn.o
global_adaptive_algorithm.mod: global_all.o

# Object file dependencies
check.o: precision_model.mod internal_types.mod
cui.o: precision_model.mod internal_types.mod check_input.mod error_handling.mod ds_routines.mod cubaturerule_general.mod global_adaptive_algorithm.mod
divide.o: precision_model.mod internal_types.mod
ds_routines.o: precision_model.mod
global_all.o: precision_model.mod internal_types.mod volume_computation.mod region_processor.mod ds_routines.mod
internal_types.o: precision_model.mod
region_processor.o: precision_model.mod internal_types.mod subdivisions.mod  cubaturerule_general.mod
rule_1.o: precision_model.mod
rule_c2.o: precision_model.mod
rule_c3.o: precision_model.mod
rule_cn.o: precision_model.mod
rule_general.o: precision_model.mod internal_types.mod quadraturerule.mod cubaturerule_c2.mod cubaturerule_c3.mod cubaturerule_cn.mod cubaturerule_t2.mod cubaturerule_t3.mod cubaturerule_tn.mod
rule_t2.o: precision_model.mod internal_types.mod
rule_t3.o: precision_model.mod
rule_tn.o: precision_model.mod
volume.o: precision_model.mod

##############################################################################
libcubpack.a: $(OBJS_CUBPACK)
	$(AR) rvu libcubpack.a $(OBJS_CUBPACK)


$(OBJS_CUBPACK):%.o:%.f90
	$(FC_PREAMBLE) $(FC_WRAPPER) $(FFLAGS) $(MPI_FLAGS) $(OPENMP_FLAGS) -c $<


$(MODS_CUBPACK):%.mod:%.o
	

distclean: clean
	/bin/rm -f *.a 


clean:
	/bin/rm -f *.o *.mod *__genmod.f90


dump_variables:
	@echo MODS_CUBPACK=$(MODS_CUBPACK)
	@echo QUALIKIZ=$(QUALIKIZ)
