# Makefile for compiling the code in ConLaw folder

OBJECTS = functions.o flux.o godunov.o
OBJECTS2 = functions.o flux.o lax.o
MODULES = functions.mod flux.mod

FFLAGS = -g

.PHONY : godunov clean

godunov : godunov.exe
	./godunov.exe

make lax : lax.exe
	./lax.exe

godunov.exe : $(MODULES) $(OBJECTS)
	gfortran $(FFLAGS) $(OBJECTS) -o godunov.exe

lax.exe : $(MODULES) $(OBJECTS2)
	gfortran $(FFLAGS) $(OBJECTS2) -o lax.exe

%.mod : %.f90
	gfortran $(FFLAGS) -c $<

%.o : %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm -rf -v *.mod *.o *.exe *.txt
