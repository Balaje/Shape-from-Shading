# $UWHPSC/codes/fortran/newton/Makefile

OBJECTS = newton.o prados_implicit.o
MODULES = newton.mod

FFLAGS = -g

.PHONY: prados_implicit clean

prados_implicit: prados_implicit.exe
	./prados_implicit.exe

prados_implicit.exe: $(MODULES) $(OBJECTS)
	gfortran -fdefault-real-8 -ffree-line-length-none $(FFLAGS) $(OBJECTS) -o prados_implicit.exe

%.o : %.f90
	gfortran -fdefault-real-8 -ffree-line-length-none $(FFLAGS) -c  $<

%.mod: %.f90
	gfortran -fdefault-real-8 -ffree-line-length-none $(FFLAGS) -c $<

clean:
	rm -f *.o *.exe *.mod op_neumannf.txt
