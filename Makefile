QLKDIR?=$(abspath .)
QUALIKIZ_SRC?=$(QLKDIR)/src
include $(QUALIKIZ_SRC)/Makefile.inc
QUALIKIZ_SRC?=$(QLKDIR)/src
QUALIKIZ_LIBSRC?=$(QLKDIR)/lib/src
CUBPACK_DIR?=$(QUALIKIZ_LIBSRC)/cubpack
GENZ_DIR?=$(QUALIKIZ_LIBSRC)/genz
NAG_DIR?=$(QUALIKIZ_LIBSRC)/nag
SLATEC_DIR?=$(QUALIKIZ_LIBSRC)/slatec
SPECFUN_DIR?=$(QUALIKIZ_LIBSRC)/specfun
FUKUSHIMA_DIR?=$(QUALIKIZ_LIBSRC)/fukushima

CUBPACK_LIB?=$(CUBPACK_DIR)/libcubpack.a
GENZ_LIB?=$(GENZ_DIR)/libgenz.a
NAG_LIB?=$(NAG_DIR)/libnag.a
SLATEC_LIB?=$(SLATEC_DIR)/libslatec.a
SPECFUN_LIB?=$(SPECFUN_DIR)/libspecfun.a
FUKUSHIMA_LIB?=$(FUKUSHIMA_DIR)/libfukushima.a

QUALIKIZ_LIBS?=$(CUBPACK_LIB) $(GENZ_LIB) $(NAG_LIB) $(SLATEC_LIB) $(SPECFUN_LIB) $(FUKUSHIMA_LIB)
QUALIKIZ_LIBS_CLEAN?=$(QUALIKIZ_LIBS:%=%.clean)

##############################################################################
QualiKiz: $(QUALIKIZ_LIBS)
	make -C $(QUALIKIZ_SRC) QuaLiKiz
	cp -f $(QUALIKIZ_SRC)/QuaLiKiz .


$(LIBNAME): $(QUALIKIZ_LIBS) $(QUALIKIZ_SRC)/Makefile.inc
	make -C $(QUALIKIZ_SRC) qlk_tci_module.mod
	ar vr $(LIBNAME) $? qlk_tci_module.o nanfilter.o

$(QUALIKIZ_SRC)/Makefile.inc:
	cp $(QUALIKIZ_SRC)/make.inc/Makefile.jetto $(QUALIKIZ_SRC)/Makefile.inc


libs: $(QUALIKIZ_LIBS)


$(QUALIKIZ_LIBS):
	make -C $(@D) $(@F)


$(QUALIKIZ_LIBS_CLEAN):
	@echo cleaning $(@D)
	-make -C $(@D) distclean


clean: $(QUALIKIZ_LIBS_CLEAN)
	make -C $(QUALIKIZ_SRC) distclean


distclean realclean: clean
	rm -f QuaLiKiz


dump_variables:
	@echo QUALIKIZ_SRC=$(QUALIKIZ_SRC)
	@echo QUALIKIZ_LIBS=$(QUALIKIZ_LIBS)
	@echo QUALIKIZ_LIBS_CLEAN=$(QUALIKIZ_LIBS_CLEAN)
	@echo QLKDIR=$(QLKDIR)
