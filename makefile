
OptFlag=-g
PetscPath=${PETSC}include
PetscBuildPath=${PETSC}arch_WSL_Build_A/include
Include=-I"${PetscPath}" -I"${PetscBuildPath}"
Module=
Libs=-lpetsc -llapack
FC=mpif90

Flags:=-cpp ${OptFlag} ${Include} ${Module}

all:hello.exe test_petsc.exe

globals.o: globals.f90
	${FC}  -c globals.f90  ${Flags}
fem_order2.o: fem_order2.f90
	${FC}  -c fem_order2.f90 ${Flags}
set.o: set.f90
	${FC}  -c set.f90 ${Flags}
para_csr.o: para_csr.f90
	${FC} -c para_csr.f90 ${Flags}

hello.exe: hello.f90 readgrid2.f90 globals.o fem_order2.o set.o para_csr.o
	${FC}  readgrid2.f90 globals.o fem_order2.o set.o para_csr.o hello.f90 -o hello.exe  ${Flags}  ${Libs}

test_petsc.exe: test_petsc.f90
	${FC}  test_petsc.f90 ${Flags} ${Libs} -o test_petsc.exe

test_set.exe: test_set.f90 set.f90
	${FC} test_set.f90 set.f90 ${Flags} -o test_set.exe

.PHONY:clean

clean:
	rm *.o *.exe
