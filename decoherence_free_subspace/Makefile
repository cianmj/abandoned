goal: dfs.x

.SUFFIXES: .f90 .o

IFOR = ifort -FR -g -w95
# -O3
EXTRALIB = -llapack

OBJECTS = dfs_main.o dfs_smatrix.o dfs_parameters.o

dfs.x : $(OBJECTS)
	$(IFOR) $(OBJECTS) $(EXTRALIB) -o dfs.x

dfs_main.o: dfs_main.f90 dfs_smatrix.o dfs_parameters.o

dfs_smatrix.o: dfs_smatrix.f90 dfs_parameters.o

dfs_parameters.o : dfs_parameters.f90

.f90.o:
	$(IFOR) -c $<

clean:
	-rm *.x *.mod *.out *.o
