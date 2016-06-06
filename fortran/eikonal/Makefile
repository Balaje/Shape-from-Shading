# Makefile for compiling the fortran code

OBJECTS = func.o god.o
MODULES = func.mod
OBJECTS2 = func.o god1.o

FFLAGS = -g

.PHONY: god clean

convex : god.exe
	./god.exe

concave : god1.exe
	./god1.exe

god.exe : $(OBJECTS) $(MODULES)
	gfortran $(FFLAGS) $(OBJECTS) -o god.exe

god1.exe : $(OBJECTS2) $(MODULES)
	gfortran $(FFLAGS) $(OBJECTS2) -o god1.exe

%.o : %.f90
	gfortran $(FFLAGS) -c $<

%.mod : %.f90
	gfortran $(FFLAGS) -c $<

clean :
	rm -f *.o *.mod *.exe *.txt
