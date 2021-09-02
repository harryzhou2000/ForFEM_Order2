
OptFlag=-g -O3
PetscPath=${PETSC_DIR}/include
PetscBuildPath=${PETSC_DIR}/${PETSC_ARCH}/include
SlepcPath=${SLEPC_DIR}/include
SlepcBuildPath=${SLEPC_DIR}/${PETSC_ARCH}/include
Include=-I"${PetscPath}" -I"${PetscBuildPath}"  -I"${SlepcPath}" -I"${SlepcBuildPath}"
Module=
Libs=-lslepc -lpetsc -lsuperlu_dist -llapack -lparmetis -lmetis -lstdc++
FC=mpiifort

Flags:=-cpp ${OptFlag} ${Include} ${Module}
FortranTargets:=globals.o set.o para_csr.o  common_utils.o elastic_constitution.o fem_order2.o

all: main_cooler.exe test_petsc.exe

# para_csr.o: para_csr.f90
# 	${FC} -c $^  -o $@ ${Flags}
# globals.o: globals.f90
# 	${FC} -c $^  -o $@ ${Flags}
# fem_order2.o: fem_order2.f90
# 	${FC} -c $^  -o $@ ${Flags}
# set.o: set.f90
# 	${FC} -c $^  -o $@ ${Flags}
# common_utils.o: common_utils.f90
# 	${FC} -c $^  -o $@ ${Flags}
# elastic_constitution.o: elastic_constitution.f90
# 	${FC} -c $^  -o $@ ${Flags}

%.o: %.f90
	${FC} -c $^  -o $@ ${Flags}


main.exe: main.f90 readgrid2.f90 ${FortranTargets}
	${FC} $^  -o $@  ${Flags}  ${Libs}
main_cooler.exe: main_cooler.f90 readgrid2.f90 ${FortranTargets}
	${FC} $^  -o $@  ${Flags}  ${Libs}
main_beam.exe: main_beam.f90 readgrid2.f90 ${FortranTargets}
	${FC} $^  -o $@  ${Flags}  ${Libs}

test_petsc.exe: test_petsc.f90
	${FC}  test_petsc.f90 ${Flags} ${Libs} -o test_petsc.exe

test_set.exe: test_set.f90 set.f90 common_utils.f90
	${FC} common_utils.f90 test_set.f90 set.f90 ${Flags} ${Libs} -o test_set.exe

.PHONY:clean

clean:
	rm -f *.o *.exe *.mod
