# edit FC and FCFLAGS as appropriate to your machine
FC = gfortran
FCFLAGS = -g -O0

bxsf2scad.x : bxsf2scad.o polyhedron.o
	$(FC) $(FCFLAGS) -o bxsf2scad.x $^

bxsf2scad.o : polyhedron.o bxsf2scad.f90
	$(FC) $(FCFLAGS) -c -o bxsf2scad.o bxsf2scad.f90

polyhedron.o : polyhedron.f90
	$(FC) $(FCFLAGS) -c -o polyhedron.o polyhedron.f90

clean:
	-rm -rf *.o *.x *.mod *.x.dSYM
