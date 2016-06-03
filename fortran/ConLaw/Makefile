# Makefile for compiling the code in ConLaw folder

OBJECTS = functions.o flux.o godunov.o
MODULES = functions.mod flux.mod

FFLAGS = -g

.PHONY : godunov clean

godunov : godunov.exe
	./godunov.exe

godunov.exe : $(MODULES) $(OBJECTS)
	gfortran $(FFLAGS) $(OBJECTS) -o godunov.exe

%.mod : %.f90
	gfortran $(FFLAGS) -c $<

%.o : %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm -rf -v *.mod *.o *.exe *.txt
