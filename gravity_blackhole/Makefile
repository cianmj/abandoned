goal: run.x

.SUFFIXES: .f90 .o
 
IFOR = gfortran
EXTRALIB = /usr/lib/liblapack.so.3gf
EXTRALIB2 = /usr/lib/libblas.so.3gf

OBJECTS = grsolve.o tnewton.o tclague.o quartic.o

run.x:	$(OBJECTS)
	$(IFOR) $(OBJECTS) $(EXTRALIB) $(EXTRALIB2) -o run.x

grsolve.o: grsolve.f90 tnewton.o tclague.o quartic.o

tnewton.o: tnewton.f90

tclague.o: tclague.f90

quartic.o: quartic.f90

.f90.o:
	$(IFOR) -c $<

clean:
	-rm *.x *.mod fort.* *.o
