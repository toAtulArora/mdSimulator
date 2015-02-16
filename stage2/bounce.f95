program bouncer
use gnuplot_fortran

implicit none
real :: dt=0.001
!N is the number of particles
!mN is the number of iterations for which a movie is made
!tN is the number of itrations (not a parameter because we want to keep it variable)
!NOTE: You must initialize Nmax with the maximum number of particles you wish to simulate the system with | else you'll get memory overflow errrors
integer(kind=4) :: N=1000
integer(kind=4), parameter :: tN=10000, mN=500,Nmax=1000000

!this is the volume
real :: volume=10
!this is the mass
real, parameter :: m = 1


real :: boxSize
real :: area
real :: time=0
!t is time
!E is energy
!P is pressure (instaneous)
!Pw is windowed pressure (averaged over optimal window size (0.1 percent)
real, dimension(tN) :: t,E,P,P2,P3,P4,Pw

integer, dimension(10) :: windowSize
!integer :: pWindowSize = 1
real :: temp1
real :: avgP
!there're six walls, two normal to each of the 3 axis
real, dimension(3,2) :: pressureWall

real, dimension(3,tN) :: qT,qDotT
real, dimension(3,Nmax) :: q
real, dimension(3,Nmax) :: qDot
integer(kind=4) :: k,i,j
!real, dimension(1) :: plX,plY,plZ
!put some number as seed


write(*,*) "the largest integer i can store is ", huge(i)
write (*,*) "initializing initial q and qDot"
call init()
write (*,*) "done"
call startPlot()

call setXrange(0.0,real(boxSize))
call setYrange(0.0,real(boxSize))
call setZrange(0.0,real(boxSize))

write (*,*) "Starting iterations.."
call iterateTheseManyTimes(tN,.true.)
write (*,*) "Done iterating :)"
windowSize(1)=evaluateAvgWindowSizeForP(Pw)

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
call plot2dSave(t( 1:(tN-windowSize(1)) ),Pw ( 1:(tN-windowSize(1)) ),"averagedP.jpg") 

!call plot2dSave(t(1:(tN-i)),Pw(1:(tN-i)),rangeYstart=0.0,rangeYend=(avgP*2.0),rangeXstart=0.0,rangeXend=(N*dt)) 

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

  function signedRand()
    real :: signedRand
    signedRand=rand()
    if(rand()>0.5) then
       signedRand=-signedRand
    end if
  end function signedRand

  !To initialize the position and velocity of the particles with random values
  subroutine init()
    !intialize boxSize etc. required for iterating
    boxSize=volume**(1.0/3.0)
    area=volume**(2.0/3.0)

    !box size is between 0 and 1
    call srand(1000)

    !initialize particles to be at some random location
    E(1)=0
    t(:)=0

    do i=1,N
       q(1,i)=rand()*boxSize
       q(2,i)=rand()*boxSize
       q(3,i)=rand()*boxSize

       qDot(1,i)=signedRand()*500 !rand()*10
       qDot(2,i)=signedRand()*500 !rand()*10
       qDot(3,i)=signedRand()*500  !rand()*10
       E(1)=E(1)+energy(qDot(:,i))


       !if(mod(i,1000)==1) then
          !write(*,'(A)',advance='no') "#"
          !write(*,'(I5.5A)',advance='no') i,", "
       !end if
    end do
  end subroutine init

  subroutine iterateTheseManyTimes(numberOfIterations,plotGraphs)
    integer(kind=4), intent(in) :: numberOfIterations
    logical, optional :: plotGraphs
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
       qT(:,k)=q(:,1)
       qDotT(:,k)=qDot(:,1)
       P(k)=sum(sum(pressureWall,dim=1),dim=1)
       if (present(plotGraphs)) then
          if(k<mN) then
             call nextPlot3d(q(1,:),q(2,:),q(3,:))
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
