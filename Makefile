#FF = ifort -xHost -Ofast
#MPIFF = mpiifort -xHost -Ofast
#LIBS =  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#OUT = a.i.out

# GNU
FF = gfortran -O2 -ffree-form
LIBS = -lblas -llapack   # on a laptop
OUT= a.out

SOURCES = det_n2.f event_mod.f rndm_mumbers.f

all:
	$(FF) $(SOURCES) fermi-hubbard.f $(LIBS) -o $(OUT)

mpi:
	$(MPIFF) $(SOURCES) fermi-hubbard-mpi.f $(LIBS) -o a.i.out.mpi

clean:
	rm a.out *mod
