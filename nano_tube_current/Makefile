goal: well.x

.SUFFIXES: .f90 .o

IFOR = ifort -FR -w95 -g
EXTRALIB = -llapack

OBJECTS = 1dwell.o subr.o potential.o param.o

well.x:	$(OBJECTS)
	$(IFOR) $(OBJECTS) $(EXTRALIB) -o well.x

1dwell.o: 1dwell.f90 subr.o potential.o param.o

subr.o: subr.f90 param.o

potential.o: potential.f90 param.o

param.o: param.f90

.f90.o:
	$(IFOR) -c $<

clean:
	-rm *.x *.mod fort.* *.o
