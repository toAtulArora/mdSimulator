!STAGE 3 modernized
program bouncer
use gnuplot_fortran
use random
use stat
implicit none
!old dt=0.001
real :: dt=0.001, radius=0.01
!radius is the radius of the hard sphere
!N is the number of particles
!mN is the number of iterations for which a movie is made
!tN is the number of itrations (not a parameter because we want to keep it variable)
!NOTE: You must initialize Nmax with the maximum number of particles you wish to simulate the system with | else you'll get memory overflow errrors
type macroscopicLike
   integer(kind=4) :: N=1000000
   !this is the volume
   !real :: volume=1000
   real:: V=1000
   !time average P and time average E
   !real :: avgP, avgE
   real :: P, E
end type macroscopicLike

type(macroscopicLike) :: macroscopic

!old tN=3000
integer(kind=4), parameter :: tN=3000, mN=1,Nmax=1000000
integer, parameter :: collisionAlgorithm = 1
!this is the mass
real, parameter :: m = 1
!this is the mass of 1 gass molecule (in Kg)
!real, parameter :: m=1.673534e10-24
!boltzman constant in joules per kelvin
real, parameter :: Kb=1.3806488e-23

integer :: tEquiv=0


real :: boxSize
real :: area
real :: time=0

!forgot what this does | probably something to do with finidng the time for averagein or some such
integer, dimension(10000) :: windowSize,windowSizeN

!integer :: pWindowSize = 1
real :: temp1,temp2

type particleLike
   real, dimension(3) :: q
   real, dimension(3) :: qDot
end type particleLike


type frameLike
   !there're six walls, two normal to each of the 3 axis
   real, dimension(3,2) :: pressureWall
   type(particleLike), dimension(Nmax) ::  particles
end type frameLike

type(frameLike) :: frame,lastFrame


type frameRecordLike
!t is time
!E is energy
!P is pressure (instaneous)
!Pw is windowed pressure (averaged over optimal window size (0.1 percent)
   real :: t,E,P,P2,P3,P4,Pw
   real, dimension(3) :: qT,qDotT
end type frameRecordLike

type (frameRecordLike), dimension(tN) :: framesRecord




! type (particle), dimension(Nmax) :: particles

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
l=32-8
!do l=1,48
!l=1
open(unit=4,file='PVminusNKTvsN')
open(unit=5,file='PavgVsN')

do tEquiv=1,50
   macroscopic%N=10**(l/8.0)
   write (*,*) "Initializing initial q and qDot for ", macroscopic%N, " particles, with tEquiv ", tEquiv
   call init(real(1.2**tEquiv))

   write (*,*) "Done!\n"

   !histOutput = hist(qDot(2,1:N),10)
   !write(*,*) "THIS IS WHAT I SENT: ",qDot(2,1:N)
   !call nextPlot2d(histOutput(1:10,1),histOutput(1:10,2))

   !write(*,*) qDot(2,1:N)

   call setXrange(0.0,real(boxSize))
   call setYrange(0.0,real(boxSize))
   call setZrange(0.0,real(boxSize))

   write (*,*) "Starting iterations (simulating ",dt*tN," physical seconds)..."
   lastFrame=frame
   call cpu_time(temp1)   
   call iterateTheseManyTimes(tN,collisionAlgorithm,1)
   call cpu_time(temp2)
   write(*,*) "Done iterating in  ", temp2- temp1, " seconds\n"
   !write (*,*) "Done iterating :) \n"
   !l=l+1
   
      
   macroscopic%P=sum(framesRecord(1000:tN)%P)/real(tN-1000)
   !volume is known
   macroscopic%E=sum(framesRecord(1000:tN)%E)/real(tN-1000)   
   !avgE=sum(E(1000:tN))/real(tN-1000)
   !temp1=avgP*volume - (N*Kb*temperature(avgE))
   
   temp1=macroscopic%P*macroscopic%V - (2/3.0)*macroscopic%E
   write(4,*) macroscopic%N,temp1,macroscopic%P,macroscopic%E,macroscopic%P*macroscopic%V,(2.0/3.0)*macroscopic%E




   !windowSize(l)=evaluateAvgWindowSizeForP(Pw,15.0)
   !windowSizeN(l)=N

   !histOutput = hist(qDot(2,1:N),10)
   !call nextPlot2d(histOutput(1:10,1),histOutput(1:10,2))

   !write(*,*) q(2,1:N)
   !E(1:tN)
   histOutput=hist(framesRecord(1:tN)%E,50)
   !write(*,*) E(1:tN)
   !call nextPlot2d(histOutput(1:50,1),histOutput(1:50,2))
   !call plot2dSave(histOutput(1:50,1),histOutput(1:50,2),"EnergyHist.jpg")
   
!,rangeXstart=0,rangeXend=6e7)
   !avgE=sum(E)/real(tN)
   !write(4,*) N,stdDevFromHist(histOutput(1:50,:)),stdDevFromHist(histOutput(1:50,:))/avgP

   !histogram of x component of qDot

end do

close(4)
close(5)
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
call plot2dSave(framesRecord(1:tN)%t,framesRecord(1:tN)%qT(3),"qT3.pdf",picFormat=1)
call plot2dSave(framesRecord(1:tN)%t,framesRecord(1:tN)%qT(2),"qT2.pdf",picFormat=1)
call plot2dSave(framesRecord(1:tN)%t,framesRecord(1:tN)%qT(1),"qT1.pdf",picFormat=1)
!  call plot2dSave(t,qT(2,:),"qT2.pdf")
!  call plot2dSave(t,qT(3,:),"qT3.pdf")
!  call plot2dSave(t,E,"energyT.pdf")
call plot2dSave(framesRecord(1:tN)%t,framesRecord(1:tN)%qDotT(3),"qDotT3.pdf",picFormat=1)
call plot2dSave(framesRecord(1:tN)%t,framesRecord(1:tN)%qDotT(2),"qDotT2.pdf",picFormat=1)
call plot2dSave(framesRecord(1:tN)%t,framesRecord(1:tN)%qDotT(1),"qDotT1.pdf",picFormat=1)
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

  !function removeElementFromArray
  
  !To initialize the position and velocity of the particles with random values
  subroutine init(parametricT)
    real, optional:: parametricT
    !intialize boxSize etc. required for iterating
    boxSize=macroscopic%V**(1.0/3.0)
    area=macroscopic%V**(2.0/3.0)
    
    !box size is between 0 and 1
    !call srand(1000)
    
    !initialize particles to be at some random location
    framesRecord(1)%E=0
    framesRecord(:)%t=0

    do i=1,macroscopic%N
       !Its the spooknfiguration
       ! q(1,i)=((real(i)/real(N))*(boxSize-(2*radius))) + radius
       ! q(2,i)=2*radius
       ! q(3,i)=2*radius

       frame%particles(i)%q(1)=rand()*(boxSize-(2*radius)) + radius
       frame%particles(i)%q(2)=rand()*(boxSize-(2*radius)) + radius
       frame%particles(i)%q(3)=rand()*(boxSize-(2*radius)) + radius

       if (present(parametricT)) then
          frame%particles(i)%qDot(1)=randomNormal()*parametricT !rand()*10
          frame%particles(i)%qDot(2)=randomNormal()*parametricT !rand()*10
          frame%particles(i)%qDot(3)=randomNormal()*parametricT  !rand()*10
       else
          !on an average, the velocity will be 0,
          !max velocity will be 500
          frame%particles(i)%qDot(1)=randomNormal()*1000 !rand()*10
          frame%particles(i)%qDot(2)=randomNormal()*1000 !rand()*10
          frame%particles(i)%qDot(3)=randomNormal()*1000  !rand()*10
       end if


       ! qDot(1,i)=100 !randomNormal()*10000 !rand()*10
       ! qDot(2,i)=50 !randomNormal()*100000 !rand()*10
       ! qDot(3,i)=10 !randomNormal()*100000  !rand()*10

       ! qDot(1,i)=signedRand()*100 !rand()*10
       ! qDot(2,i)=signedRand()*100 !rand()*10
       ! qDot(3,i)=signedRand()*100  !rand()*10
       framesRecord(1)%E=framesRecord(1)%E+energy(frame%particles(i)%qDot(:))


       !if(mod(i,1000)==1) then
          !write(*,'(A)',advance='no') "#"
          !write(*,'(I5.5A)',advance='no') i,", "
       !end if
    end do
  end subroutine init

  subroutine nMatPrettyPrint(nMat)
    !real :: nMatPrettyPrint
    integer, dimension(:,:,:,:)::nMat
    integer :: sizePrec,sizeDensity,i,j
    sizePrec=size(nMat(:,1,1,1))
    sizeDensity=size(nMat(1,1,1,:))
    write(*,*) "nMat=["
    do i=1,sizePrec
       write(*,*) nMat(:,i,1,1)
    end do
    write(*,*) "]"    
  end subroutine nMatPrettyPrint





  subroutine iterateTheseManyTimes(numberOfIterations,collisionSpheres,plotGraphs)
    integer(kind=4), intent(in) :: numberOfIterations
    integer, optional :: plotGraphs,collisionSpheres
    !integer :: i,ii,j
    integer(kind=4) :: ii
    !this is for efficient collision detection
    integer(kind=4), dimension(3) :: qDisc
    real, dimension(3) :: qVector,qDotVector,unitqVector,qDotNewI,qDotNewII
    logical :: collision
    !the a,b,c are for solving a quadratic! heheh
    real :: deltaT,a,b,c
    !this is for neighbour searching
    integer, parameter :: nPrec=20,density=10
    integer :: allocateStatus
    integer, dimension(:,:,:,:), allocatable :: nMat
    !integer(kind=4), dimension(0:nPrec,0:nPrec,0:nPrec,0:density)::nMat=0
    integer(kind=4), dimension (0:density) :: nMatCurr=0
    integer(kind=4) :: densityInBox,ni
    !call nMatPrettyPrint(nMat)
    call initProgress()
    !write(*,'(A)') "\b[##-------]"
    allocate(nMat(0:nPrec,0:nPrec,0:nPrec,0:density),stat=allocateStatus)
    if(allocateStatus /= 0) stop "Not enough memory :("
    do k=1,numberOfIterations
       framesRecord(k)%t=time
       time=time+dt
       framesRecord(k)%E=0
       framesRecord(k)%P=0
       frame%pressureWall(:,:)=0
       !collision with the walls
       !initialize neigbhour matrix to zero
       nMat=0
       !write(*,*) "Gla2"
       do i=1,macroscopic%N
          !update particle position
          ! q(:,i)=q(:,i)+(qDot(:,i)*dt)
          frame%particles(i)%q(:)=frame%particles(i)%q(:)+(frame%particles(i)%qDot(:)*dt)

          !check for wall collisions and evaluate pressure on each wall on the fly
          do j=1,3
             !if (q(j,i) >= 1.0) then
             !if (q(j,i) >= boxSize) then
             if (frame%particles(i)%q(j)>=boxSize-radius) then
                frame%pressureWall(j,1)=frame%pressureWall(j,1)+ ( (2*m*frame%particles(i)%qDot(j))/ (dt*area) )
                frame%particles(i)%qDot(j)=-frame%particles(i)%qDot(i)
                frame%particles(i)%q(j)=frame%particles(i)%q(j)-2*(frame%particles(i)%q(j)-(boxSize-radius))
             !else if (q(j,i) <= 0.0) then
             else if (frame%particles(i)%q(j)<=radius) then
                frame%pressureWall(j,2)=frame%pressureWall(j,2)-( (2*m*frame%particles(i)%qDot(j)) / (dt*area) )
                frame%particles(i)%qDot(j)=-frame%particles(i)%qDot(j)
                frame%particles(i)%q(j) = (2*radius)-frame%particles(i)%q(j) 
             end if                          
          end do
          !write(*,"(AI4.4)",advance="no") "\rWallCollisionsDone:",i
          !evaluate total energy
          framesRecord(k)%E=framesRecord(k)%E+energy(frame%particles(i)%qDot(:))


          !check for collisions as hard spheres
          if(present(collisionSpheres)) then
             !write(*,*) "Will check for collisions using"
             if(collisionSpheres==1) then
                !write(*,*) "Fast algorithm"
                !First discritize the particle's location
                qDisc=int((frame%particles(i)%q(:)/boxSize)*nPrec)
                !extract the relavent box from the neighbour matrix
                !This will contain the number of particles in side the box 
                !(that aren't colliding among themselves by construction)
                !if(qDisc(1).ge.0 .and. qDisc(1).le.nPrec .and. qDisc(2).ge.0 .and. qDisc(2).le.nPrec .and. qDisc(2).ge.0 .and. qDisc(2).le.nPrec) then
                   nMatCurr=nMat(qDisc(1),qDisc(2),qDisc(3),:)
                !else
                   !write (*,*) qDisc(:)
                   !stop "ERROR!"
                !end if

                !IF the box was empty, coolio :)
                if(nMatCurr(0)==0) then
                   !tell it that there's a particle
                   nMatCurr(0)=1
                   !and tell it which particle
                   nMatCurr(1)=i
                   !Now if the box is NOT empty, then
                else
                   !find out how many particles are already there
                   densityInBox=nMatCurr(0)
                   !ERROR, but let it run
                   if(densityInBox>density) then
                      call errorHandle("Too many particles in one grid")
                      densityInBox=density
                   end if
                   
                   !this is to find out if the ith particle collides with 
                   !the ones already there in the box
                   collision=.false.
                   !we loop over all particles inside the box
                   do ni=1,densityInBox
                      !test for collisions between the i and nMatCurr(ni) particles

                      !COLLISION TEST AND RESPONSE ALGO START:
                      ii=nMatCurr(ni)
                      qVector=frame%particles(i)%q(:) - frame%particles(ii)%q(:)                      
                      qDotVector=frame%particles(i)%qDot(:) - frame%particles(ii)%qDot(:)

                      !if the distance between the centres is less than 2r
                      !then COLLISION HAPPENED!
                      !TODO: Figure what to do when the idstance is exactly zero! :(
                      if (lenSquare(qVector)<4*(radius*radius)) then
                         if (lenSquare(qVector)>0) then
                            unitqVector=qVector/length(qVector)
                            !delta T is the time elapsed since the actual collision
                            !Courtesy Prashansa
                            a=lenSquare(qDotVector)
                            b=-2*(sum(qVector*qDotVector))
                            c=lenSquare(qVector) -4*radius*radius
                            if (a>0) then
                               deltaT=(-b + sqrt((b*b) - (4*a*c)))/(2*a)
                               !Update position: use old velocity to get back to the point of collision
                               frame%particles(i)%q(:)=frame%particles(i)%q(:) - (frame%particles(i)%qDot(:)*deltaT)
                               frame%particles(ii)%q(:)=frame%particles(ii)%q(:) - (frame%particles(ii)%qDot(:)*deltaT)

                               !update velocity

                               qDotNewI=frame%particles(i)%qDot(:)-sum(unitqVector*frame%particles(i)%qDot(:))*unitqVector + sum(unitqVector*frame%particles(ii)%qDot(:))*unitqVector
                               !qDotNewI=qDotNewI+sum(unitqVector*frame%particles(ii)%qDot(:))*unitqVector
                               !had to split because else it won't compile :(
                               qDotNewII=frame%particles(ii)%qDot(:)-sum(unitqVector*frame%particles(ii)%qDot(:))*unitqVector + sum(unitqVector*frame%particles(i)%qDot(:))*unitqVector 
                               !qDotNewII=qDotNewII+sum(unitqVector*frame%particles(i)%qDot(:))*unitqVector

                               frame%particles(i)%qDot(:)=qDotNewI
                               frame%particles(ii)%qDot(:)=qDotNewII

                               !Update position: use corrected velocity to update to where you should've been, 
                               !had the collision been detected when it happened! Damange control..hehe
                               frame%particles(i)%q(:)=frame%particles(i)%q(:) + (frame%particles(i)%qDot(:)*deltaT)
                               frame%particles(ii)%q(:)=frame%particles(ii)%q(:) + (frame%particles(ii)%qDot(:)*deltaT)

                            else

                               !This shouldn't happen! It means that two particles are moving in the same direction and they're closer than 2r
                               !So we fix it by hand..
                               frame%particles(ii)%q(:)=frame%particles(i)%q(:)+ 2*radius*unitqVector
                               !leave the volicities unchanged
                            end if
                         end if
                         !COLLISION TEST AND RESPONSE ALGO: ENDED

                         !OBVIOUSLY THERE WAS A COLLISION, SO DO THE NEEDFUL FOR THE NEIGHBOUR SEARCH ALGO
                         !make the number of particles in the box decrement
                         !since by assumption, there wont be any 3 body collisions
                         nMatCurr(0)=nMatCurr(0)-1
                         !remove ni th element from the nMat matrix                   
                         nMatCurr=(/ nMatCurr(0:(ni-1)),nMatCurr((ni+1):density) /)
                         !make collision happened true for use later
                         collision=.true.
                         !exit the loop (by assumption, there can't be anymore collisions)
                         exit
                      end if
                   end do
                   !meaning no collision was found
                   if(collision .eqv. .false.) then
                      !increment the number of particles in the box that're not self colliding
                      nMatCurr(0)=nMatCurr(0)+1
                      !add the address of the current particle in the list of particles in the box
                      nMatCurr(nMatCurr(0))=i
                   end if
                end if

                !update the nMat matrix
                nMat(qDisc(1),qDisc(2),qDisc(3),:)=nMatCurr
                !will print the matrix once the loop is done
                ! if (i==N .and. k<mN) then
                !    call nMatPrettyPrint(nMat)
                ! end if
             end if
          end if
       end do

       !collision among spheres
       if (present(collisionSpheres)) then
          if (collisionSpheres==2) then
          !write(*,*) "Collision test using algorithm 2"
          do i=1,macroscopic%N
             do ii=i+1,macroscopic%N
             !if(i .ne. ii) then

                collision=.false.
                !weak test for collision
                do j=1,3
                   if (abs(frame%particles(i)%q(j)-frame%particles(ii)%q(j))<2*radius) then
                      collision=.true.
                      exit
                   end if
                end do
                !if weak test says collision, test more carefully
                if (collision) then
                   qVector=frame%particles(i)%q(:)-frame%particles(ii)%q(:)
                   qDotVector=frame%particles(i)%qDot(:) - frame%particles(ii)%qDot(:)

                   !if the distance between the centres is less than 2r
                   !then COLLISION HAPPENED!               
                   if (lenSquare(qVector)<4*(radius*radius)) then

                      !delta T is the time elapsed since the actual collision
                      !Courtesy Prashansa
                      a=lenSquare(qDotVector)
                      b=-2*(sum(qVector*qDotVector))
                      c=lenSquare(qVector) -4*radius*radius
                      deltaT=(-b + sqrt((b*b) - (4*a*c)))/(2*a)

                      !Update position: use old velocity to get back to the point of collision
                      frame%particles(i)%q(:)=frame%particles(i)%q(:) - (frame%particles(i)%qDot(:)*deltaT)
                      frame%particles(ii)%q(:)=frame%particles(ii)%q(:) - (frame%particles(ii)%qDot(:)*deltaT)

                      !update velocity
                      unitqVector=qVector/length(qVector)
                      qDotNewI=frame%particles(i)%qDot(:)-sum(unitqVector*frame%particles(i)%qDot(:))*unitqVector+sum(unitqVector*frame%particles(ii)%qDot(:))*unitqVector
                      ! qDotNewI=qDotNewI+sum(unitqVector*frame%particles(ii)%qDot(:))*unitqVector
                      qDotNewII=frame%particles(ii)%qDot(:)-sum(unitqVector*frame%particles(ii)%qDot(:))*unitqVector+sum(unitqVector*frame%particles(i)%qDot(:))*unitqVector
                      ! qDotNewII=qDotNewII+sum(unitqVector*frame%particles(i)%qDot(:))*unitqVector

                      frame%particles(i)%qDot(:)=qDotNewI
                      frame%particles(ii)%qDot(:)=qDotNewII

                      !Update position: use corrected velocity to update to where you should've been, had the collision been detected when it happened! Damange control..hehe
                      frame%particles(i)%q(:)=frame%particles(i)%q(:) + (frame%particles(i)%qDot(:)*deltaT)
                      frame%particles(ii)%q(:)=frame%particles(ii)%q(:) + (frame%particles(ii)%qDot(:)*deltaT)
                      
                   end if

                end if
             !end if
             end do
          end do
          end if
       end if

       !THis is to save the location and position of particle 1
       !at different times
       framesRecord(k)%qT(:)=frame%particles(1)%q(:)
       framesRecord(k)%qDotT(:)=frame%particles(1)%qDot(:)
       !average pressure over all walls (or equivalently, could've put the surface area as 6*area of square)
       framesRecord(k)%P=sum(sum(frame%pressureWall,dim=1),dim=1)/6.0
       if (present(plotGraphs)) then
          if(k<mN) then
             if(plotGraphs==1 .or. plotGraphs==3) then                
                call nextPlot3d(frame%particles(1:macroscopic%N)%q(1),frame%particles(1:macroscopic%N)%q(2),frame%particles(1:macroscopic%N)%q(3))!q(2,1:N),q(3,1:N))
             end if
             if(plotGraphs==2 .or. plotGraphs==3) then
                
                histOutput = hist((frame%particles(1:macroscopic%N)%qDot(1)**2) + (frame%particles(1:macroscopic%N)%qDot(2)**2) + (frame%particles(1:macroscopic%N)%qDot(3)**2),20)
                ! hist((qDot(2,1:N)*qDot(2,1:N))+ (qDot(1,1:N)*qDot(1,1:N)) + (qDot(3,1:N)*qDot(3,1:N)),20)
                call nextPlot2d(histOutput(1:20,1),histOutput(1:20,2))
             end if
          end if
       end if
       !write(*,'(f5.5A)',advance="no") real(k)/real(numberOfIterations),"\r"
       call incProgress(k,numberOfIterations)
    end do
    deallocate(nMat)
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

  subroutine errorHandle(errorMessage)    
    character(len=*) :: errorMessage
    integer :: inp
    write(*,*) "Error Occured:"
    write(*,*) errorMessage
    write(*,*) "Input 0 to ignore and continue, anything else to exit"
    ! read(*,*) inp
    ! if (inp .ne. 0) then
       stop "User chose to exit"
    ! end if
  end subroutine errorHandle

  function avgPwithTheseManyIterations(window)
    integer(kind=4),intent(in) :: window
    real, dimension(tN) :: avgPwithTheseManyIterations
    integer :: h,k,l,m,n
    real :: tempP
    do h=1,tN-window
       tempP=0
       do k=h,(h+window)
          tempP=tempP + framesRecord(k)%P
       end do
       tempP=tempP/real(window)
       avgPwithTheseManyIterations(h)=tempP
    end do
  end function avgPwithTheseManyIterations

  function errorPercentInPn(Pn)
    real, intent(in), dimension(tN) :: Pn
    real :: stddev,errorPercentInPn
    stddev = sqrt( sum((/Pn-macroscopic%P/)*(/Pn-macroscopic%P/),dim=1)   / (tN-1))
    errorPercentInPn = (stddev/macroscopic%P)*100
    
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

    macroscopic%P=sum(framesRecord(1:tN)%P)/real(tN)
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
