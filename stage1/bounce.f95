program bouncer
use gnuplot_fortran
implicit none
real :: dt=0.001,t=0
real, dimension(10) :: x,y,z
real, dimension(10) :: vx,vy,vz
integer :: N=size(x),k,i
!real, dimension(1) :: plX,plY,plZ
!put some number as seed

!box size is between 0 and 1
call srand(1000)

!initialize particles to be at some random location
do i=1,N
   x(i)=rand()
   y(i)=rand()
   z(i)=rand()
   
   vx(i)=rand()*10
   vy(i)=rand()*10
   vz(i)=rand()*10
end do

!make x array
!dx = (x_end - x_start)/n        
!x(0:n) = [(i*dx, i=0,n)]        


call startPlot()
do k=1,100
   do i=1,N
      t=t+dt
      x(i)=x(i)+(vx(i)*dt)
      y(i)=y(i)+(vy(i)*dt)
      z(i)=z(i)+(vz(i)*dt)
   end do
   call nextPlot3d(x,y,z)
end do

!do i=1,100
   !make y array
   !y=sin(x) / (i*x*(1.0/100.0)+1)
   !call nextPlot3d(x,y,x)
!end do

call endPlot()


end program bouncer
