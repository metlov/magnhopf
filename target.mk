.SUFFIXES:

ifndef _ARCH
_ARCH := $(shell uname -m)
export _ARCH
endif

OBJDIR := _bin_$(_ARCH)

MAKETARGET = [ -d $(OBJDIR) ] && \
             $(MAKE) -I $(OBJDIR) --no-print-directory -C $(OBJDIR) \
                 -f $(CURDIR)/Makefile \
                 SRCDIR=$(CURDIR) $(MAKECMDGOALS)

.PHONY: $(OBJDIR)
$(OBJDIR):
	+@[ -d $@ ] || mkdir -p $@
	+@$(MAKETARGET)

Makefile : ;
%.mk :: ;

% :: $(OBJDIR) ; :

.PHONY: clean
clean:
	+@$(MAKETARGET) || true
	rm -rf $(OBJDIR)
