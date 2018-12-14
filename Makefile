#!/bin/make -f
QUALIKIZ_SRC?=src
include $(QUALIKIZ_SRC)/Makefile.inc
QUALIKIZ_SRC?=$(QUALIKIZ)/src
QUALIKIZ_LIBSRC?=$(QUALIKIZ)/lib/src
SLATEC_DIR?=$(QUALIKIZ_LIBSRC)/slatec
SPECFUN_DIR?=$(QUALIKIZ_LIBSRC)/specfun
FUKUSHIMA_DIR?=$(QUALIKIZ_LIBSRC)/fukushima

SLATEC_LIB?=$(SLATEC_DIR)/libslatec.a
SPECFUN_LIB?=$(SPECFUN_DIR)/libspecfun.a
FUKUSHIMA_LIB?=$(FUKUSHIMA_DIR)/libfukushima.a

QUALIKIZ_LIBS?=$(SLATEC_LIB) $(SPECFUN_LIB) $(FUKUSHIMA_LIB)
QUALIKIZ_LIBS_CLEAN?=$(QUALIKIZ_LIBS:%=%.clean)

##############################################################################
QualiKiz: $(QUALIKIZ_LIBS)
	make -C $(QUALIKIZ_SRC) QuaLiKiz
	cp -f $(QUALIKIZ_SRC)/QuaLiKiz .

integration_driver: $(QUALIKIZ_LIBS)
	make -C $(QUALIKIZ_SRC) integration_driver
	cp -f $(QUALIKIZ_SRC)/integration_driver .

libs: $(QUALIKIZ_LIBS)


$(QUALIKIZ_LIBS):
	make -C $(@D) $(@F)


$(QUALIKIZ_LIBS_CLEAN):
	@echo cleaning $(@D)
	-make -C $(@D) distclean


clean: $(QUALIKIZ_LIBS_CLEAN)
	make -C $(QUALIKIZ_SRC) distclean


distclean: clean
	rm -f QuaLiKiz


dump_variables:
	@echo QUALIKIZ_SRC=$(QUALIKIZ_SRC)
	@echo QUALIKIZ_LIBS=$(QUALIKIZ_LIBS)
	@echo QUALIKIZ_LIBS_CLEAN=$(QUALIKIZ_LIBS_CLEAN)
	@echo QUALIKIZ=$(QUALIKIZ)
