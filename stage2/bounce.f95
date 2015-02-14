program bouncer
use gnuplot_fortran

implicit none
real :: dt=0.001
!N is the number of particles
!tN is the number of itrations
!mN is the number of iterations for which a movie is made
integer, parameter :: N=10000, tN=5000, mN=500
!this is the volume
real, parameter :: volume=10
!this is the mass
real, parameter :: m = 1

real, parameter :: boxSize=volume**(1.0/3.0)
real, parameter :: area = volume**(2.0/3.0)
real :: time=0
real, dimension(tN) :: t,E,P,P2,P3,P4,PnTemp
integer :: pWindowSize = 1
real :: temp1
real :: avgP
real, dimension(3,2) :: pressureWall

real, dimension(3,tN) :: qT,qDotT
real, dimension(3,N) :: q
real, dimension(3,N) :: qDot
integer :: k,i,j
!real, dimension(1) :: plX,plY,plZ
!put some number as seed




call init()

call startPlot()

call setXrange(0.0,real(boxSize))
call setYrange(0.0,real(boxSize))
call setZrange(0.0,real(boxSize))

call iterateTheseManyTimes(tN,.true.)

avgP=sum(P)/tN
do i=1,100
   PnTemp=avgPwithTheseManyIterations(i)
   temp1=errorPercentInPn(PnTemp)
   !write(*,*) "Error in P",i,": ",temp1
   if (temp1<0.1) then
      write(*,*) "Window size is ",i
      write(*,*) "Time over which average is required is ",real(i)*dt, " and the error then is ",temp1
      pWindowSize=i
      exit
   end if
   
   !call plot2dSave(t(1:size(PnTemp)),PnTemp,"presAv")
end do

!PnTemp

! P2 = avgPwithTheseManyIterations(2)
! write(*,*) "error in P2:", errorPercentInPn(P2)
! P3 = avgPwithTheseManyIterations(3)
! write(*,*) "error in P3:", errorPercentInPn(P4)
! P4 = avgPwithTheseManyIterations(10)
! write(*,*) "error in P4:", errorPercentInPn(P4)

call plot2dSave(t,qT(1,:),"qT1.jpg")
call plot2dSave(t,qT(2,:),"qT2.jpg")
call plot2dSave(t,qT(3,:),"qT3.jpg")
call plot2dSave(t,E,"energyT.jpg")
call plot2dSave(t,qDotT(1,:),"qDotT1.jpg")
call plot2dSave(t,qDotT(2,:),"qDotT2.jpg")
call plot2dSave(t,qDotT(3,:),"qDotT3.jpg")
call plot2dSave(t,P,"pressureT.jpg")
call plot2dSave(t(1:size(PnTemp)),PnTemp,"averagedP.jpg")
! call plot2dSave(t,P2,"pressure2T.jpg",rangeYstart=100000.0,rangeYend=550000.0)
! call plot2dSave(t,P3,"pressure3T.jpg",rangeYstart=100000.0,rangeYend=550000.0)
! call plot2dSave(t,P4,"pressure4T.jpg",rangeYstart=100000.0,rangeYend=550000.0)


call endPlot()


contains
  function energy(x)
    real, dimension(3), intent(in) :: x
    real :: energy
    energy = ( (x(1)*x(1)) + (x(2)*x(2)) + (x(3)*x(3)) ) * (m/2.0)
  end function energy

  !To initialize the position and velocity of the particles with random values
  subroutine init()
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
  end subroutine init

  subroutine iterateTheseManyTimes(numberOfIterations,plotGraphs)
    integer, intent(in) :: numberOfIterations
    logical, optional :: plotGraphs
    do k=1,numberOfIterations
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
       end do
       qT(:,k)=q(:,1)
       qDotT(:,k)=qDot(:,1)
       P(k)=sum(sum(pressureWall,dim=1),dim=1)
       if (present(plotGraphs)) then
          if(k<mN) then
             call nextPlot3d(q(1,:),q(2,:),q(3,:))
          end if
       end if
    end do
  end subroutine iterateTheseManyTimes

  function avgPwithTheseManyIterations(window)
    integer,intent(in) :: window
    real, dimension(tN) :: avgPwithTheseManyIterations
    integer :: h,k,l,m,n
    real :: tempP
    do h=1,tN-window
       tempP=0
       do k=h,(h+window)
          tempP=tempP + P(k)
       end do
       tempP=tempP/real(window)
       avgPwithTheseManyIterations(h)=tempP
    end do
  end function avgPwithTheseManyIterations

  function errorPercentInPn(Pn)
    real, intent(in), dimension(tN) :: Pn
    real :: stddev,errorPercentInPn
    stddev = sqrt( sum((/Pn-avgP/)*(/Pn-avgP/),dim=1)  ) / tN
    errorPercentInPn = (stddev/avgP)*100
    
  end function errorPercentInPn
    

end program bouncer


!make x array
!dx = (x_end - x_start)/n        
!x(0:n) = [(i*dx, i=0,n)]        
