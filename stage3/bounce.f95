program bouncer
use gnuplot_fortran
use random
use stat
implicit none
real :: dt=0.001, radius=0.5
!radius is the radius of the hard sphere
!N is the number of particles
!mN is the number of iterations for which a movie is made
!tN is the number of itrations (not a parameter because we want to keep it variable)
!NOTE: You must initialize Nmax with the maximum number of particles you wish to simulate the system with | else you'll get memory overflow errrors
integer(kind=4) :: N=5
integer(kind=4), parameter :: tN=1000, mN=500,Nmax=1000000

!this is the volume
real :: volume=1000
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

integer, dimension(10000) :: windowSize,windowSizeN
!integer :: pWindowSize = 1
real :: temp1
real :: avgP
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


!radius = 4.0

call startPlot()


!write(*,*) "the largest integer i can store is ", huge(i)

!l=0
!N=1000
!open(unit=4,file='energyHist')

!do l=1,48
!l=1

   !N=10**(l/8.0)
   write (*,*) "Initializing initial q and qDot for ", N, " particles."
   call init()
   write (*,*) "Done!\n"

   !histOutput = hist(qDot(2,1:N),10)
   !write(*,*) "THIS IS WHAT I SENT: ",qDot(2,1:N)
   !call nextPlot2d(histOutput(1:10,1),histOutput(1:10,2))

   !write(*,*) qDot(2,1:N)

   call setXrange(0.0,real(boxSize))
   call setYrange(0.0,real(boxSize))
   call setZrange(0.0,real(boxSize))

   write (*,*) "Starting iterations.."

   call iterateTheseManyTimes(tN,1,1)
   write (*,*) "Done iterating :) \n"
   !l=l+1
   
   !windowSize(l)=evaluateAvgWindowSizeForP(Pw,15.0)
   !windowSizeN(l)=N

   !histOutput = hist(qDot(2,1:N),10)
   !call nextPlot2d(histOutput(1:10,1),histOutput(1:10,2))

   !write(*,*) q(2,1:N)

   histOutput=hist(E(1:tN),50)
   !write(*,*) E(1:tN)
   !call nextPlot2d(histOutput(1:50,1),histOutput(1:50,2))
   !call plot2dSave(histOutput(1:50,1),histOutput(1:50,2),"EnergyHist.jpg")
   
!,rangeXstart=0,rangeXend=6e7)
   !avgE=sum(E)/real(tN)
   !write(4,*) N,stdDevFromHist(histOutput(1:50,:)),stdDevFromHist(histOutput(1:50,:))/avgP

   !histogram of x component of qDot

!end do

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
call plot2dSave(t(1:tN),qT(3,1:tN),"qT3.pdf",picFormat=1)
call plot2dSave(t(1:tN),qT(2,1:tN),"qT2.pdf",picFormat=1)
call plot2dSave(t(1:tN),qT(1,1:tN),"qT1.pdf",picFormat=1)
!  call plot2dSave(t,qT(2,:),"qT2.pdf")
!  call plot2dSave(t,qT(3,:),"qT3.pdf")
!  call plot2dSave(t,E,"energyT.pdf")
call plot2dSave(t(1:tN),qDotT(3,1:tN),"qDotT3.pdf",picFormat=1)
call plot2dSave(t(1:tN),qDotT(2,1:tN),"qDotT2.pdf",picFormat=1)
call plot2dSave(t(1:tN),qDotT(1,1:tN),"qDotT1.pdf",picFormat=1)
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

  function signedRand()
    real :: signedRand
    signedRand=rand()
    if(rand()>0.5) then
       signedRand=-signedRand
    end if
  end function signedRand

  function lenSquare(vector)
    real, dimension(3):: vector
    real :: lenSquare
    lenSquare=vector(1)*vector(1) + vector(2)*vector(2) + vector(3)*vector(3)
  end function lenSquare

  function length(vector)
    real :: length
    real, dimension(3) :: vector
    length=sqrt(lenSquare(vector))
  end function length

  
  !To initialize the position and velocity of the particles with random values
  subroutine init()
    !intialize boxSize etc. required for iterating
    boxSize=volume**(1.0/3.0)
    area=volume**(2.0/3.0)
    
    !box size is between 0 and 1
    !call srand(1000)
    
    !initialize particles to be at some random location
    E(1)=0
    t(:)=0

    do i=1,N
       q(1,i)=((real(i)/real(N))*(boxSize-(2*radius))) + radius
       q(2,i)=2*radius
       q(3,i)=2*radius
       !q(1,i)=rand()*(boxSize-(2*radius)) + radius
       !q(2,i)=rand()*(boxSize-(2*radius)) + radius
       !q(3,i)=rand()*(boxSize-(2*radius)) + radius

       qDot(1,i)=randomNormal()*100000 !rand()*10
       qDot(2,i)=randomNormal()*100000 !rand()*10
       qDot(3,i)=randomNormal()*100000  !rand()*10

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

  subroutine iterateTheseManyTimes(numberOfIterations,collisionSpheres,plotGraphs)
    integer(kind=4), intent(in) :: numberOfIterations
    integer, optional :: plotGraphs,collisionSpheres
    !integer :: i,ii,j
    integer :: ii
    real, dimension(3) :: qVector,qDotVector,unitqVector,qDotNewI,qDotNewII
    logical :: collision
    !the a,b,c are for solving a quadratic! heheh
    real :: deltaT,a,b,c
    call initProgress()
    !write(*,'(A)') "\b[##-------]"
    do k=1,numberOfIterations
       t(k)=time
       time=time+dt
       E(k)=0
       P(k)=0
       pressureWall(:,:)=0
       !collision with the walls
       do i=1,N
          q(:,i)=q(:,i)+(qDot(:,i)*dt)
          do j=1,3
             !if (q(j,i) >= 1.0) then
             !if (q(j,i) >= boxSize) then
             if (q(j,i)>=boxSize-radius) then
                pressureWall(j,1)=pressureWall(j,1)+ ( (2*m*qDot(j,i))/ (dt*area) )
                qDot(j,i)=-qDot(j,i)
                q(j,i)=q(j,i)-2*(q(j,i)-(boxSize-radius))
             !else if (q(j,i) <= 0.0) then
             else if (q(j,i)<=radius) then
                pressureWall(j,2)=pressureWall(j,2)-( (2*m*qDot(j,i)) / (dt*area) )
                qDot(j,i)=-qDot(j,i)
                q(j,i) = (2*radius)-q(j,i) 
             end if
          end do
          E(k)=E(k)+energy(qDot(:,i))
       end do

       !collision among spheres
       if (present(collisionSpheres)) then
          if (collisionSpheres==1) then
          do i=1,N
             do ii=i+1,N
             !if(i .ne. ii) then

                collision=.false.
                !weak test for collision
                do j=1,3
                   if (abs(q(j,i)-q(j,ii))<2*radius) then
                      collision=.true.
                      exit
                   end if
                end do
                !if weak test says collision, test more carefully
                qVector=q(:,i)-q(:,ii)
                qDotVector=qDot(:,i) - qDot(:,ii)
                if (collision) then
                   !if the distance between the centres is less than 2r
                   !then COLLISION HAPPENED!               
                   if (lenSquare(qVector)<4*(radius*radius)) then
                      !write(*,*) "There was a collision:" 

                      !delta T is (2r-r')/v_rel
                      !deltaT=((2*radius) - length(qVector))/length(qDotVector)
                      !deltaT=(2.0*radius)/length(qDotVector)

                      !delta T is the time elapsed since the actual collision
                      !Courtesy Prashansa
                      a=lenSquare(qDotVector)
                      b=-2*(sum(qVector*qDotVector))
                      c=lenSquare(qVector) -4*radius*radius
                      deltaT=(-b + sqrt((b*b) - (4*a*c)))/(2*a)

                      !write(*,*) "V1", qDot(:,i) ," and V2", qDot(:,ii)
                      !write(*,*) "It had actually heppened (:P) ", deltaT, " unit times ago"
                      !nonsense done
                      ! deltaT=((2*radius) - length(qVector))/length(qDot(:,i))
                      ! q(:,i)=q(:,i) - sum(unitqVector(:)*qDot(:,i))*unitqVector(:)*deltaT
                      ! deltaT=((2*radius) - length(qVector))/length(qDot(:,ii))                      
                      ! q(:,ii)=q(:,ii) - sum(unitqVector(:)*qDot(:,ii))*unitqVector(:)*deltaT
                      !end nonsense

                      !Update position: use old velocity to get back to the point of collision
                      q(:,i)=q(:,i) - (qDot(:,i)*deltaT)
                      q(:,ii)=q(:,ii) - (qDot(:,ii)*deltaT)

                      !update velocity

                      !write(*,*) "Velocity before", qDot(:,i)
                      unitqVector=qVector/length(qVector)
                      qDotNewI=qDot(:,i)-sum(unitqVector*qDot(:,i))*unitqVector+sum(unitqVector*qDot(:,ii))*unitqVector
                      qDotNewII=qDot(:,ii)-sum(unitqVector*qDot(:,ii))*unitqVector+sum(unitqVector*qDot(:,i))*unitqVector
                      ! qDot(:,i)=qDot(:,i) - 2*sum(unitqVector*qDot(:,i))*unitqVector 
                      ! qDot(:,ii)=qDot(:,ii) - 2*sum(unitqVector*qDot(:,ii))*unitqVector

                      qDot(:,i)=qDotNewI
                      qDot(:,ii)=qDotNewII

                      !write(*,*) "Velocity after: ", qDot(:,i)
                      !qDotVector

                      !Update position: use corrected velocity to update to where you should've been, had the collision been detected when it happened! Damange control..hehe
                      q(:,i)=q(:,i) + (qDot(:,i)*deltaT)
                      q(:,ii)=q(:,ii) + (qDot(:,ii)*deltaT)
                      
                   end if

                end if
             !end if
             end do
          end do
          end if
       end if

       !THis is to save the location and position of particle 1
       !at different times
       qT(:,k)=q(:,1)
       qDotT(:,k)=qDot(:,1)
       P(k)=sum(sum(pressureWall,dim=1),dim=1)
       if (present(plotGraphs)) then
          if(k<mN) then
             if(plotGraphs==1 .or. plotGraphs==3) then

                call nextPlot3d(q(1,1:N),q(2,1:N),q(3,1:N))
             end if
             if(plotGraphs==2 .or. plotGraphs==3) then
                histOutput = hist((qDot(2,1:N)*qDot(2,1:N))+ (qDot(1,1:N)*qDot(1,1:N)) + (qDot(3,1:N)*qDot(3,1:N)),20)
                call nextPlot2d(histOutput(1:20,1),histOutput(1:20,2))
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
