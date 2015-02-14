program bouncer
use gnuplot_fortran

implicit none
real :: dt=0.001
!N is the number of particles
!tN is the number of itrations
!mN is the number of iterations for which a movie is made
integer, parameter :: N=1000, tN=5000, mN=500
!this is the volume
real, parameter :: volume=10
!this is the mass
real, parameter :: m = 1

real, parameter :: boxSize=volume**(1.0/3.0)
real, parameter :: area = volume**(2.0/3.0)
real :: time=0
real, dimension(tN) :: t,E,P
real, dimension(3,2) :: pressureWall

real, dimension(3,tN) :: qT,qDotT
real, dimension(3,N) :: q
real, dimension(3,N) :: qDot
integer :: k,i,j
!real, dimension(1) :: plX,plY,plZ
!put some number as seed


!box size is between 0 and 1
call srand(1000)

!initialize particles to be at some random location
E(1)=0
do i=1,N
   t(i)=0
   q(1,i)=rand()*boxSize
   q(2,i)=rand()*boxSize
   q(3,i)=rand()*boxSize
   
   qDot(1,i)=10 !rand()*10
   qDot(2,i)=20 !rand()*10
   qDot(3,i)=30 !rand()*10
   E(1)=E(1)+energy(qDot(:,i))
end do

!make x array
!dx = (x_end - x_start)/n        
!x(0:n) = [(i*dx, i=0,n)]        


call startPlot()

call setXrange(0.0,real(boxSize))
call setYrange(0.0,real(boxSize))
call setZrange(0.0,real(boxSize))


do k=1,tN
   t(k)=time
   time=time+dt
   E(k)=0
   P(k)=0
   pressureWall(:,:)=0
   do i=1,N
      q(:,i)=q(:,i)+(qDot(:,i)*dt)
      do j=1,3
         !if (q(j,i) >= 1.0) then
         if (q(j,i) >= boxSize) then
            pressureWall(j,1)=pressureWall(j,1)+ ( (2*m*qDot(j,i))/ (dt*area) )
            qDot(j,i)=-qDot(j,i)
            q(j,i)=q(j,i)-2*(q(j,i)-boxSize)
         else if (q(j,i) <= 0.0) then
            pressureWall(j,2)=pressureWall(j,2)-( (2*m*qDot(j,i)) / (dt*area) )
            qDot(j,i)=-qDot(j,i)
            q(j,i) = -q(j,i) 
         end if
      end do
      E(k)=E(k)+energy(qDot(:,i))
      !y(i)=y(i)+(vy(i)*dt)
      !z(i)=z(i)+(vz(i)*dt)
      ! if (x(i) > 1.0) then
      !    vx(i)=-vx(i)
      !    x(i)=x(i)- 2*( 1-x(i))
      ! end if
      ! if (y(i) > 1.0) then
      !    vy(i)=-vy(i)
      !    y(i)=y(i)- 2*( 1-y(i))
      ! end if
      ! if (z(i) > 1.0) then
      !    vz(i)=-vz(i)
      !    z(i)=z(i)- 2*( 1-z(i))
      ! end if
   end do
   qT(:,k)=q(:,1)
   qDotT(:,k)=qDot(:,1)
   P(k)=sum(sum(pressureWall,dim=1),dim=1)
   
   if(k<mN) then
      call nextPlot3d(q(1,:),q(2,:),q(3,:))
   end if
end do
call plot2dSave(t,qT(1,:),"qT1.jpg")
call plot2dSave(t,qT(2,:),"qT2.jpg")
call plot2dSave(t,qT(3,:),"qT3.jpg")
call plot2dSave(t,E,"energyT.jpg")
call plot2dSave(t,qDotT(1,:),"qDotT1.jpg")
call plot2dSave(t,qDotT(2,:),"qDotT2.jpg")
call plot2dSave(t,qDotT(3,:),"qDotT3.jpg")
call plot2dSave(t,P,"pressureT.jpg")

!do i=1,100
   !make y array
   !y=sin(x) / (i*x*(1.0/100.0)+1)
   !call nextPlot3d(x,y,x)
!end do

call endPlot()


contains
  function energy(x)
    real, dimension(3), intent(in) :: x
    real :: energy
    energy = ( (x(1)*x(1)) + (x(2)*x(2)) + (x(3)*x(3)) ) * (m/2.0)
  end function energy

end program bouncer
