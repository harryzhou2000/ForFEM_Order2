
OptFlag=-g -O0
PetscPath=${PETSC}include
PetscBuildPath=${PETSC}arch_WSL_Build_A/include
Include=-I"${PetscPath}" -I"${PetscBuildPath}"
Module=
Libs=-lpetsc -llapack
FC=mpif90

Flags:=-cpp ${OptFlag} ${Include} ${Module}
FortranTargets:=globals.o fem_order2.o set.o para_csr.o common_utils.o elastic_constitution.o

all:main.exe test_petsc.exe

globals.o: globals.f90
	${FC} -c $^  -o $@ ${Flags}
fem_order2.o: fem_order2.f90
	${FC} -c $^  -o $@ ${Flags}
set.o: set.f90
	${FC} -c $^  -o $@ ${Flags}
para_csr.o: para_csr.f90
	${FC} -c $^  -o $@ ${Flags}
common_utils.o: common_utils.f90
	${FC} -c $^  -o $@ ${Flags}
elastic_constitution.o: elastic_constitution.f90
	${FC} -c $^  -o $@ ${Flags}



main.exe: main.f90 readgrid2.f90 ${FortranTargets}
	${FC} $^  -o $@  ${Flags}  ${Libs}

test_petsc.exe: test_petsc.f90
	${FC}  test_petsc.f90 ${Flags} ${Libs} -o test_petsc.exe

test_set.exe: test_set.f90 set.f90
	${FC} test_set.f90 set.f90 ${Flags} -o test_set.exe

.PHONY:clean

clean:
	rm -f *.o *.exe
