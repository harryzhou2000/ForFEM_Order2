
hello.exe: hello.f90 readgrid2.f90 fem_order2.f90
	gfortran readgrid2.f90 globals.f90 fem_order2.f90 hello.f90 -o hello.exe -g