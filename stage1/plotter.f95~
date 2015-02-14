program plotter
use gnuplot_fortran
implicit none
!character(40) :: fileName
!character*30 fileName
integer, parameter :: n = 100
real, dimension(0:n) :: x,y
real :: x_start = 0.0, x_end = 20, dx
integer :: i

!make x array
dx = (x_end - x_start)/n
x(0:n) = [(i*dx, i=0,n)]
call startPlot()

do i=1,1000
   !make y array
   y=sin(x) / (i*x*(1.0/100.0)+1)
   call nextPlot(x,y)
end do

call endPlot()
! generate data for plot
!i=1
!filename="Hello"
!write(*,*) filename
!write (fileName,'(a,i4.4,a)') 'file',i,'.jpg'
!write (*,*) fileName
!call plot2d(x,y,fileName)

!call system ("gnuplot command")
! call system ("gnuplot")
! call system ("set terminal jpeg")
! call system ("set output 'file2.jpg'")
! call system ("plot data.dat")
! call system ("exit")

end program plotter
