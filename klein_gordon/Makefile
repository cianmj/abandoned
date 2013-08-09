goal: well.x

.SUFFIXES: .f90 .o

IFOR = ifort -O3
#EXTRALIB = -llapack

OBJECTS = 1dwell.o kg_fd.o potential.o param.o

well.x:	$(OBJECTS)
	$(IFOR) $(OBJECTS) -o well.x
#	$(IFOR) $(OBJECTS) $(EXTRALIB) -o well.x

1dwell.o: 1dwell.f90 kg_fd.o potential.o param.o

kg_fd.o: kg_fd.f90 param.o

potential.o: potential.f90 param.o

param.o: param.f90

.f90.o:
	$(IFOR) -c $<

clean:
	-rm *.x *.mod *.o
