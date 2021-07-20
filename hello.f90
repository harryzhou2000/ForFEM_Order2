

program addNumbers
   use globals
   use fem_order2
   implicit none
   real :: a, b, result 
   integer :: c;
   integer, dimension(5) :: someInts
   character(80) :: outfile
   character(10) :: outtitle = "goodstart"
   real(8) :: D(3,3) = reshape((/1, 2, 3, 5 ,5, 6, 1, 8, 9/),(/3,3/))

   FILEINP = "./mark2_external.neu"
   outfile = "./out2.plt"
   call readgfile

   call output_plt_mesh(outfile, outtitle)


   a = -12.5
   b = 15.0
   result = a + b
   c = a
   someInts(1) = 1;
   someInts(2:4) = (/1,2,3/)
   call initializeLib
   call getVolumes
   call output_plt_scalar("./out2_data1.plt", "goodstart",cell_volumes,"cell_volume")

   
   print*,D
   print*, matmul(directInverse3x3(D),D)

   print *, 'The total is ', result , ' cis ', c
   
end program addNumbers

