OPT = -ffpe-trap=invalid,zero,overflow

all: aerosol_box.exe

aerosol_box.exe: aerosol_box.f90 aerosol_mod.o
	gfortran $(OPT) $^ -o $@

aerosol_mod.o: aerosol_mod.f90
	gfortran $(OPT) -c $<

clean:
	rm *.exe
	rm *.o
	rm *.mod
