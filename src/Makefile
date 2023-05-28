FC = gfortran
FCFLAGS = -g -fbounds-check -fbacktrace \
          -ffpe-trap=zero,invalid,overflow,underflow
FCLEGACY = -std=legacy
FC90FLAGS = -Wall -Wunused-parameter -Wdo-subscript -fcheck=all -std=f2008 \
            -Wno-unused-dummy-argument
CC = gcc
CFLAGS = -g

PROGRAMS = TEST PHB_TRACE HOPFSHELL

# the main computing engine
HOPFLIB_OBJS= d1mach.o xerror.o pnpoly.o \
              vdqng.o dchfdv.o dpchfd.o mus1.o mus2.o mus3.o zeroin.o \
              debug.o tracer.o quadratures.o energy.o profile.o \
              skyrmions.o groundstate.o \
              skyrme_state.o hopf_state.o

TEST_OBJS=${HOPFLIB_OBJS} test_state.o test.o
TEST_DEPS=${TEST_OBJS:.o=.d}

PHB_TRACE_OBJS=${HOPFLIB_OBJS} test_state.o phb_trace.o
PHB_TRACE_DEPS=${PHB_TRACE_OBJS:.o=.d}

HOPFSHELL_OBJS=FCreadline.o jsu_readline.o hopfshell.o ${HOPFLIB_OBJS}
HOPFSHELL_DEPS=${HOPFSHELL_OBJS:.o=.d}
HOPFSHELL_LIBS=-lreadline

OBJS = $(TEST_OBJS)
DEPENDS = $(TEST_DEPS) $(PHB_TRACE_DEPS) $(HOPFSHELL_DEPS)

all: $(PROGRAMS)

TEST: $(TEST_OBJS) $(TEST_DEPS)
	$(FC) $(FCFLAGS) $(FC90FLAGS) -o $@ $(TEST_OBJS) $(LDFLAGS)

PHB_TRACE: $(PHB_TRACE_OBJS) $(PHB_TRACE_DEPS)
	$(FC) $(FCFLAGS) $(FC90FLAGS) -o $@ $(PHB_TRACE_OBJS) $(LDFLAGS)

HOPFSHELL: $(HOPFSHELL_OBJS) $(HOPFSHELL_DEPS)
	$(FC) $(FCFLAGS) $(FC90FLAGS) -o $@ $(HOPFSHELL_OBJS) \
            $(HOPFSHELL_LIBS) $(LDFLAGS)

# note the include file dependency for initial build (before .d files are
# generated)
groundstate.mod groundstate.o: skyrmion_ground.inc

%.o %d: %.f
	$(FC) $(FCFLAGS) $(FCLEGACY) -cpp -MD -c $<
	@touch $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o %.d %.mod: %.f90
	$(FC) $(FCFLAGS) $(FC90FLAGS) -cpp -MD -c $<
	@touch $@

%.o %.d %.mod: %.F90
	$(FC) $(FCFLAGS) $(FC90FLAGS) -MD -c $<
	@touch $@

%.d: %.f
	$(FC) $(FCFLAGS) $(FCLEGACY) -MD -cpp -M $< > $@

%.d: %.c
	$(CC) $(CFLAGS) -M $< -o $@

skyrmion_ground.inc: skyrmion_ground.dat dat2inc.awk
	DAT2INC=`echo $^ | cut -f2 -d' '`; \
        LINES=`awk -E $${DAT2INC} $< | wc -l`; \
        (echo "double precision :: data_q($${LINES}) = [    &"; \
         awk -f $${DAT2INC} col=1 $< ; \
         echo "]" ; \
         echo "double precision :: data_h($${LINES}) = [    &"; \
         awk -f $${DAT2INC} col=2 $< ; \
         echo "]" ; \
         ) > $(SRCDIR)/src/skyrmion_ground.inc

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f $(OBJS) $(DEPENDS)
	rm -f $(SRCDIR)/src/skyrmion_ground.inc

veryclean: clean
	rm -f *~ $(PROGRAMS)

-include $(DEPENDS)
