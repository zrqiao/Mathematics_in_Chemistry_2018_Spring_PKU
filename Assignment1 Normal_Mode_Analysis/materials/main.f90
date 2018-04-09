program main
   implicit none
   real*8 :: x(9), v, dv(9), ddv(9,9)
   logical :: fex
   character*30 :: buff, fin

   call get_command_argument(1, fin)
   if(fin=='') then
      write(*,*) 'Use the config file as the 1st command argument.'
      stop
   end if

   inquire(file=fin,exist=fex)
   if(.not.fex) then
      write(*,*) 'File "', trim(adjustl(fin)), '" not found.'
      stop
   end if

   open(10,file=fin)
   read(10,*) buff, x(1:3)
   read(10,*) buff, x(4:6)
   read(10,*) buff, x(7:9)
   close(10)

   call h2opot(v,dv,ddv,x,2)

   write(*,*) 'Energy (outputted to "ener.dat"):'
   write(*,*) v
   open(11,file='ener.dat')
   write(11,*) v
   close(11)

   write(*,*) 'Gradient (outputted to "grad.dat"):'
   write(*,'(9F)') dv
   open(12,file='grad.dat')
   write(12,'(1F)') dv
   close(12)

   write(*,*) 'Hessian (outputted to "hess.dat"):'
   write(*,'(9F)') ddv
   open(13,file='hess.dat')
   write(13,'(9F)') ddv
   close(13)
end program
