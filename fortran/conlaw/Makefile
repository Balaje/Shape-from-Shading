# Makefile for compiling the fortran code

OBJECTS = fun.o god.o
OBJECTS2 = fun.o lax.o
MODULES = fun.mod

FFLAGS = -g

.PHONY: god clean

godunov : god.exe
	./god.exe

lax : lax.exe
	./lax.exe

god.exe : $(OBJECTS) $(MODULES)
	gfortran $(FFLAGS) $(OBJECTS) -o god.exe

lax.exe : $(OBJECTS2) $(MODULES)
	gfortran $(FFLAGS) $(OBJECTS2) -o lax.exe

%.o : %.f90
	gfortran $(FFLAGS) -c $<

%.mod : %.f90
	gfortran $(FFLAGS) -c $<

clean :
	rm -f *.o *.mod *.exe
