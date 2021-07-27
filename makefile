
OptFlag=-g
PetscPath=${PETSC}include
PetscBuildPath=${PETSC}arch_WSL_Build_A/include
Include=-I"${PetscPath}" -I"${PetscBuildPath}"
Module=
Libs=-lpetsc -llappack
FC=mpif90

Flags:=-cpp ${OptFlag} ${Include} ${Module}

all:hello.exe

globals.o: globals.f90
	${FC}  -c globals.f90  ${Flags}

fem_order2.o: fem_order2.f90
	${FC}  -c fem_order2.f90 ${Flags}

hello.exe: hello.f90 readgrid2.f90 globals.o fem_order2.o
	${FC}  readgrid2.f90 globals.o fem_order2.o hello.f90 -o hello.exe  ${Flags}  ${Libs}


.PHONY:clean

clean:
	rm *.o *.exe
