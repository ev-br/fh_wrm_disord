FF = gfortran -O2
LIBS = -lblas -llapack

all:
	$(FF) det_n2.f event_mod.f rndm_mumbers.f fermi-hubbard.f $(LIBS)

clean:
	rm a.out *mod
