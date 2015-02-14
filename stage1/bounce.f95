program bouncer
use gnuplot_fortran
implicit none


integer, parameter :: n = 100
real, dimension(0:n) :: x,y
real :: x_start = 0.0, x_end = 20, dx
integer :: i

!make x array
dx = (x_end - x_start)/n
x(0:n) = [(i*dx, i=0,n)]
call startPlot()

do i=1,100
   !make y array
   y=sin(x) / (i*x*(1.0/100.0)+1)
   call nextPlot3d(x,y,x)
end do

call endPlot()


end program bouncer
