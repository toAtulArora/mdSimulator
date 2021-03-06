program bouncer
use gnuplot_fortran
use random
use stat
implicit none
real :: dt=0.001
!N is the number of particles
!mN is the number of iterations for which a movie is made
!tN is the number of itrations (not a parameter because we want to keep it variable)
!NOTE: You must initialize Nmax with the maximum number of particles you wish to simulate the system with | else you'll get memory overflow errrors
integer(kind=4) :: N=1000
integer(kind=4), parameter :: tN=100, mN=500,Nmax=1000000

!this is the volume (in m^3)
real :: volume=10
!this is the mass of 1 gass molecule (in Kg)
!real, parameter :: m=1.673534e10-24
real, parameter :: m = 1
integer :: tEquiv=0


!boltzman constant in joules per kelvin
real, parameter :: Kb=1.3806488e-23

real :: boxSize
real :: area
real :: time=0
!t is time (in seconds)
!E is energy (in Joules)
!P is pressure (instaneous) (N per m^2 = Pascals)
!Pw is windowed pressure (averaged over optimal window size (0.1 percent)
real, dimension(tN) :: t,E,P,P2,P3,P4,Pw

integer, dimension(10000) :: windowSize,windowSizeN
!integer :: pWindowSize = 1
real :: temp1
real :: avgP
real :: avgE
!there're six walls, two normal to each of the 3 axis
real, dimension(3,2) :: pressureWall

real, dimension(3,tN) :: qT,qDotT
real, dimension(3,Nmax) :: q
real, dimension(3,Nmax) :: qDot
integer(kind=4) :: k,i,j,l

real, dimension(tN,2)::histOutput

!real, dimension(1) :: plX,plY,plZ
!put some number as seed

!call seed_random_number(1)

call startPlot()


!write(*,*) "the largest integer i can store is ", huge(i)

!l=0
!N=1000
open(unit=4,file='PVminusNKTvsN')
open(unit=5,file='PavgVsN')

!do l=1,48
!l=32-8
l=32-4
!tEquiv is a paremeter that changes the temperature

do tEquiv=1,50
   N=10**(l/8.0)
   write (*,*) "Initializing initial q and qDot for ", N, " particles."
   call init(real(1.2**tEquiv))
   write (*,*) "Done!\n"

   call setXrange(0.0,real(boxSize))
   call setYrange(0.0,real(boxSize))
   call setZrange(0.0,real(boxSize))

   write (*,*) "Starting iterations.."
   call iterateTheseManyTimes(tN)
   write (*,*) "Done iterating :) \n"
   
   avgP=sum(P(1:tN))/real(tN)
   !volume is known
   !time averaged E, not E per particle!
   avgE=sum(E(1:tN))/real(tN)
   !temp1=avgP*volume - (N*Kb*temperature(avgE))
   temp1=avgP*volume - (2/3.0)*avgE
   write(4,*) N,temp1,avgP,avgE,avgP*volume,(2.0/3.0)*avgE

   ! temp1=avgP - (2/3.0)*avgE*N
   ! write(4,*) N,temp1,avgP,avgE,avgP,(2.0/3.0)*avgE*N

   !PLOTTING ERROR IN PRESSURE WITH NUMBER OF PARTICLES
   histOutput=hist(P(1:tN),50)
   call nextPlot2d(histOutput(1:50,1),histOutput(1:50,2))
   avgP=sum(P(1:tN))/real(tN)
   write(5,*) N,stdDevFromHist(histOutput(1:50,:)),stdDevFromHist(histOutput(1:50,:))/avgP

   !histogram of x component of qDot

end do

!end do

close(5)
close(4)

!l=30
!call plot2dSave(log(1.0*windowSizeN(1:l)),dt*windowSize(1:l),"windowSize vs N")

!PnTemp

! P2 = avgPwithTheseManyIterations(2)
! write(*,*) "error in P2:", errorPercentInPn(P2)
! P3 = avgPwithTheseManyIterations(3)
! write(*,*) "error in P3:", errorPercentInPn(P4)
! P4 = avgPwithTheseManyIterations(10)
! write(*,*) "error in P4:", errorPercentInPn(P4)



!NOT MAKING GRAPHS FOR NOW...
! write (*,*) "Generating graphs.."
!  call plot2dSave(t,qT(1,:),"qT1.pdf")
!  call plot2dSave(t,qT(2,:),"qT2.pdf")
!  call plot2dSave(t,qT(3,:),"qT3.pdf")
!  call plot2dSave(t,E,"energyT.pdf")
!  call plot2dSave(t,qDotT(1,:),"qDotT1.pdf")
!  call plot2dSave(t,qDotT(2,:),"qDotT2.pdf")
!  call plot2dSave(t,qDotT(3,:),"qDotT3.pdf")
! call plot2dSave(t,P,"pressureT.pdf")
! call plot2dSave(t( 1:(tN-windowSize(1)) ),Pw ( 1:(tN-windowSize(1)) ),"averagedP.pdf") 
! write (*,*) "Done \n"




!call plot2dSave(t(1:(tN-i)),Pw(1:(tN-i)),rangeYstart=0.0,rangeYend=(avgP*2.0),rangeXstart=0.0,rangeXend=(N*dt)) 

! call plot2dSave(t,P2,"pressure2T.jpg",rangeYstart=100000.0,rangeYend=550000.0)
! call plot2dSave(t,P3,"pressure3T.jpg",rangeYstart=100000.0,rangeYend=550000.0)
! call plot2dSave(t,P4,"pressure4T.jpg",rangeYstart=100000.0,rangeYend=550000.0)

write (*,*) "Finalizing visualizaiton..."
call endPlot()


write (*,*) "Done \n"

write (*,*) "---------\n All done \n"

contains
  function energy(x)
    real, dimension(3), intent(in) :: x
    real :: energy
    energy = ( (x(1)*x(1)) + (x(2)*x(2)) + (x(3)*x(3)) ) * (m/2.0)
  end function energy
  
  function temperature(totalEnergy)
    real :: temperature, totalEnergy
    temperature=((2/3.0)*(totalEnergy))/N*Kb
  end function temperature

  function signedRand()
    real :: signedRand
    signedRand=rand()
    if(rand()>0.5) then
       signedRand=-signedRand
    end if
  end function signedRand

  !To initialize the position and velocity of the particles with random values
  subroutine init(parametricT)
    real, optional:: parametricT
    !intialize boxSize etc. required for iterating
    boxSize=volume**(1.0/3.0)
    area=volume**(2.0/3.0)
    
    !box size is between 0 and 1
    !call srand(1000)
    
    !initialize particles to be at some random location
    E(1)=0
    t(:)=0

    do i=1,N
       q(1,i)=rand()*boxSize
       q(2,i)=rand()*boxSize
       q(3,i)=rand()*boxSize

       if (present(parametricT)) then
          qDot(1,i)=randomNormal()*parametricT !rand()*10
          qDot(2,i)=randomNormal()*parametricT !rand()*10
          qDot(3,i)=randomNormal()*parametricT  !rand()*10
       else
          !on an average, the velocity will be 0,
          !max velocity will be 500
          qDot(1,i)=randomNormal()*1000 !rand()*10
          qDot(2,i)=randomNormal()*1000 !rand()*10
          qDot(3,i)=randomNormal()*1000  !rand()*10
       end if
       ! qDot(1,i)=signedRand()*100 !rand()*10
       ! qDot(2,i)=signedRand()*100 !rand()*10
       ! qDot(3,i)=signedRand()*100  !rand()*10
       E(1)=E(1)+energy(qDot(:,i))


       !if(mod(i,1000)==1) then
          !write(*,'(A)',advance='no') "#"
          !write(*,'(I5.5A)',advance='no') i,", "
       !end if
    end do
  end subroutine init

  subroutine iterateTheseManyTimes(numberOfIterations,plotGraphs)
    integer(kind=4), intent(in) :: numberOfIterations
    integer, optional :: plotGraphs
    call initProgress()
    !write(*,'(A)') "\b[##-------]"
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
       !THis is to save the location and position of particle 1
       !at different times
       qT(:,k)=q(:,1)
       qDotT(:,k)=qDot(:,1)
       P(k)=sum(sum(pressureWall,dim=1),dim=1)/6.0
       if (present(plotGraphs)) then
          if(k<mN) then
             if(plotGraphs==1) then

                call nextPlot3d(q(1,:),q(2,:),q(3,:))
             else if(plotGraphs==2) then
                histOutput = hist(qDot(2,1:N),10)
                call nextPlot2d(histOutput(1:10,1),histOutput(1:10,2))
             end if
          end if
       end if
       !write(*,'(f5.5A)',advance="no") real(k)/real(numberOfIterations),"\r"
       call incProgress(k,numberOfIterations)
    end do
    call endProgress()
  end subroutine iterateTheseManyTimes

  subroutine initProgress()
    write (*,'(A)',advance="no") " [----------]\r ["
  end subroutine initProgress

  subroutine endProgress()
    write (*,*) "\r [##########]"
  end subroutine endProgress

  subroutine incProgress(currentVal,finalVal)
    integer(kind=4),intent(in) :: currentVal,finalVal
    if(mod(currentVal,int(finalVal/10.0)) ==1 ) then
       write(*,"(A)",advance="no") "#"
    end if
  end subroutine incProgress

  function avgPwithTheseManyIterations(window)
    integer(kind=4),intent(in) :: window
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
    stddev = sqrt( sum((/Pn-avgP/)*(/Pn-avgP/),dim=1)   / (tN-1))
    errorPercentInPn = (stddev/avgP)*100
    
  end function errorPercentInPn

  !PnTemp is populated with the windowed pressure
  function evaluateAvgWindowSizeForP(PnTemp,maxPercentErrorGiven)    
    real, intent(out), dimension(tN) :: PnTemp
    real, intent(in), optional :: maxPercentErrorGiven
    real :: maxPercentError=0.1
    integer::evaluateAvgWindowSizeForP
    if (present(maxPercentErrorGiven)) then
       maxPercentError=maxPercentErrorGiven
    end if

    avgP=sum(P)/real(size(P))
    do i=1,tN/10,(tN/1000+1)
       PnTemp=avgPwithTheseManyIterations(i)
       temp1=errorPercentInPn(PnTemp(1:(tN-i)) )
       write(*,*) "Error in P",i,": ",temp1
       
       if (temp1<maxPercentError) then
          write(*,*) "Window size is ",i
          write(*,*) "Time over which average is required is ",real(i)*dt, " and the error then is ",temp1
          evaluateAvgWindowSizeForP=i
          !call plot2dSave(t(1:(tN-i)),PnTemp(1:(tN-i)),rangeYstart=0.0,rangeYend=(avgP*2.0),rangeXstart=0.0,rangeXend=(tN*dt)) 
          exit
       end if
       !if(mod(i,10)==1) then          
           !call plot2dSave(t(1:(tN-i)),PnTemp(1:(tN-i)),rangeYstart=0.0,rangeYend=(avgP*2.0),rangeXstart=0.0,rangeXend=(N*dt)) 
       !end if
       !call plot2dSave(t(1:size(PnTemp)),PnTemp,"presAv")
    end do
    evaluateAvgWindowSizeForP=i
  end function evaluateAvgWindowSizeForP

end program bouncer


!make x array
!dx = (x_end - x_start)/n        
!x(0:n) = [(i*dx, i=0,n)]        
