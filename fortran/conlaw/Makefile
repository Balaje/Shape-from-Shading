# Makefile for compiling the fortran code

OBJECTS = fun.o god.o
MODULES = fun.mod

FFLAGS = -g

.PHONY: god clean

godunov : god.exe
	./god.exe

god.exe : $(OBJECTS) $(MODULES)
	gfortran $(FFLAGS) $(OBJECTS) -o god.exe

%.o : %.f90
	gfortran $(FFLAGS) -c $<

%.mod : %.f90
	gfortran $(FFLAGS) -c $<

clean :
	rm -f *.o *.mod *.exe
