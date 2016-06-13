#!/bin/make -f
SRC?=src
include $(SRC)/Makefile.inc
SRC?=$(QUALIKIZ)/src
LIBSRC?=$(QUALIKIZ)/lib/src
CUBPACKDIR?=$(LIBSRC)/cubpack
GENZDIR?=$(LIBSRC)/genz
NAGDIR?=$(LIBSRC)/nag
SLATECDIR?=$(LIBSRC)/slatec
SPECFUNDIR?=$(LIBSRC)/specfun
LIBS?=$(CUBPACKDIR)/libcubpack.a $(GENZDIR)/libgenz.a $(NAGDIR)/libnag.a $(SLATECDIR)/libslatec.a $(SPECFUNDIR)/libspecfun.a
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
