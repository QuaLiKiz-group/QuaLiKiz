#!/bin/make -f
SRC?=src
include $(SRC)/Makefile.inc
SRC?=$(QUALIKIZ)/src
LIBSRC?=$(QUALIKIZ)/lib/src
CUBPACK?=$(LIBSRC)/cubpack
GENZ?=$(LIBSRC)/genz
NAG?=$(LIBSRC)/nag
SLATEC?=$(LIBSRC)/slatec
SPECFUN?=$(LIBSRC)/specfun
LIBS?=$(CUBPACK)/libcubpack.a $(GENZ)/libgenz.a $(NAG)/libnag.a $(SLATEC)/libslatec.a $(SPECFUN)/libspecfun.a
LIBS_CLEAN?=$(LIBS:%=%.clean)

##############################################################################
QualiKiz: $(LIBS)
	make -C $(SRC) QuaLiKiz
	cp -f $(SRC)/QuaLiKiz .


libs: $(LIBS)


$(LIBS):
	make -C $(@D) $(@F)


$(LIBS_CLEAN):
	make -C $(@D) distclean


clean: $(LIBS_CLEAN)
	make -C $(SRC) distclean


distclean: clean
	rm -f QuaLiKiz


dump_variables:
	@echo LIBS=$(LIBS)
	@echo LIBS_CLEAN=$(LIBS_CLEAN)
	@echo QUALIKIZ=$(QUALIKIZ)
	@echo SRC=$(SRC)
