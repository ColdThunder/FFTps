########################################################################
#                                Makefile
#
EXEC = 3dps.out

SRCS = 3dps.f90

MODSRCS = ./head.f90 ./myomp.f90 ./p2grid.f90

MODS = $(MODSRCS:.f90=.mod)

OBJS = $(SRCS:.f90=.o) $(MODSRCS:.f90=.o)

INCL = -I./

LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

OPTS = -fPIC -shared-intel -mcmodel=large -qopenmp

FC   = ifort

$(EXEC):$(SRCS) $(MODS)
	$(FC) $(OPTS) $(SRCS) $(MODSRCS) $(INCL) $(LIBS) -o $(EXEC) 

$(MODS):$(MODSRCS)
	$(FC) -c $(OPTS) $(MODSRCS) $(INCL) $(LIBS)

.PHONY:clean
clean:
	-rm -rf $(EXEC) $(OBJS) $(MODS)
########################################################################

