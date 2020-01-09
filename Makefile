FF = gfortran -O2
MPIFF = mpif90 -O2
LIBS = -lopenblas
#LIBS = -lblas -llapack   # on a laptop

SOURCES = det_n2.f event_mod.f rndm_mumbers.f

all:
	$(FF) $(SOURCES) fermi-hubbard.f $(LIBS)

mpi:
	$(MPIFF) $(SOURCES) fermi-hubbard-mpi.f $(LIBS) -o a.out.mpi

clean:
	rm a.out *mod
