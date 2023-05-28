ifeq (,$(filter _%,$(notdir $(CURDIR))))
include target.mk
else
VPATH = $(SRCDIR)/src:$(SRCDIR)/third_party/SLATEC:$(SRCDIR)/third_party/PCHIP:\
        $(SRCDIR)/third_party/MUS $(SRCDIR)/third_party/NETLIB \
        $(SRCDIR)/third_party/JSU_READLINE $(SRCDIR)/third_party/PNPOLY \
        $(SRCDIR)/third_party/QUADPACK
include $(SRCDIR)/src/Makefile
endif
