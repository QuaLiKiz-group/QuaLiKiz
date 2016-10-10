#!/bin/make -f
SRC?=src
include $(SRC)/Makefile.inc
SRC?=$(QUALIKIZ)/src
LIBSRC?=$(QUALIKIZ)/lib/src
CUBPACK_DIR?=$(LIBSRC)/cubpack
GENZ_DIR?=$(LIBSRC)/genz
NAG_DIR?=$(LIBSRC)/nag
SLATEC_DIR?=$(LIBSRC)/slatec
SPECFUN_DIR?=$(LIBSRC)/specfun
FUKUSHIMA_DIR?=$(LIBSRC)/fukushima

CUBPACK_LIB?=$(CUBPACK_DIR)/libcubpack.a
GENZ_LIB?=$(GENZ_DIR)/libgenz.a
NAG_LIB?=$(NAG_DIR)/libnag.a
SLATEC_LIB?=$(SLATEC_DIR)/libslatec.a
SPECFUN_LIB?=$(SPECFUN_DIR)/libspecfun.a
FUKUSHIMA_LIB?=$(FUKUSHIMA_DIR)/libfukushima.a

LIBS?=$(CUBPACK_LIB) $(GENZ_LIB) $(NAG_LIB) $(SLATEC_LIB) $(SPECFUN_LIB) $(FUKUSHIMA_LIB)
LIBS_CLEAN?=$(LIBS:%=%.clean)

##############################################################################
QualiKiz: $(LIBS)
	make -C $(SRC) QuaLiKiz
	cp -f $(SRC)/QuaLiKiz .


libs: $(LIBS)


$(LIBS):
	make -C $(@D) $(@F)


$(LIBS_CLEAN):
	@echo cleaning $(@D)
	-make -C $(@D) distclean


clean: $(LIBS_CLEAN)
	make -C $(SRC) distclean


distclean: clean
	rm -f QuaLiKiz


dump_variables:
	@echo SRC=$(SRC)
	@echo LIB_SRC=$(LIBSRC)
	@echo LIBS=$(LIBS)
	@echo LIBS_CLEAN=$(LIBS_CLEAN)
	@echo QUALIKIZ=$(QUALIKIZ)
