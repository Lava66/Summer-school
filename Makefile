F90    = gfortran
FFLAGS = -O0 -g -ffpe-trap=invalid,zero,overflow -ffree-line-length-none

OBJS = main.o chemistry_mod.o opkdmain.o opkda1.o opkda2.o
EXE  = main

all: $(EXE)

main: main.f90 chemistry_mod.o opkdmain.o opkda1.o opkda2.o
	gfortran $(OPT) main.f90 chemistry_mod.o opkdmain.o opkda1.o opkda2.o -o main

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

%.o: %.f
	$(F90) $(FFLAGS) -std=legacy -w -c $< -o $@

clean:
	@rm -vf *.o *.mod main
