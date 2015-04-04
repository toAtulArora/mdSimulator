!STAGE 4

program bouncer
use gnuplot_fortran
use random
use stat
implicit none
!old dt=0.001
real,parameter :: maxVelocity=50
real,parameter :: radius=0.05,dt=0.0001 !radius/(10*maxVelocity) !0.0001
!radius is the radius of the hard sphere
!N is the number of particles
!mN is the number of iterations for which a movie is made (OBSOLETE)
!tN is the number of itrations (not a parameter because we want to keep it variable)
!NOTE: You must initialize Nmax with the maximum number of particles you wish to simulate the system with | else you'll get memory overflow errrors
type macroscopicLike
   integer(kind=4) :: N=1000
   !this is the volume
   !real :: volume=1000
   real:: V=1000
   !time average P and time average E
   !real :: avgP, avgE
   real :: P, E
end type macroscopicLike

type(macroscopicLike) :: macroscopic

!old tN=3000
real, parameter :: simulationDuration=5.0
integer(kind=4), parameter :: tN=simulationDuration/dt, Nmax=1000000,mN=100
!render duration in seconds
real, parameter :: renderDuration=5.0

integer, parameter :: collisionAlgorithm = 1
!this is the mass
real, parameter :: m = 1
real, parameter :: fps = 25.0
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

type potentialParametersLike
   real :: epsilon=100,rm=radius !0.01
   real :: A=4,B=6,X=1,Y=0
   !real :: A=12,B=6,X=1,Y=2
   !This is for the potential
   ! u = epsilon (X*(rm/r)**A - Y*(rm/r)**B)
   !d is the distance 2*radius
   !epsilon is the depth of the wander wal well
end type potentialParametersLike

type(potentialParametersLike) :: potentialParameter

type particleLike
   real, dimension(3) :: q, qDot, force, q2,q3,q4,k1,m1,k2,m2,k3,m3,k4,m4
   integer :: next
   integer :: collided
!0 = not collided
!1 = collided once and so on
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
!l=32-32
!do l=1,48
!l=1
open(unit=4,file='PVminusNKTvsN')
open(unit=5,file='PavgVsN')
! open(unit=7,file='outputForce')

! temp1 =0.04
! !write(*,*) "starting"
! do while(temp1<10.0)
!    !write(*,*) "writing"
!    write(7,*) temp1,force( (/ 0.0,0.0,0.0 /), (/ temp1, temp1,0.0/),potentialParameter)
!    temp1=temp1+0.01
! end do
! close(7)



do tEquiv=28,28

   !macroscopic%N=10**(l/8.0)
   write (*,*) "Initializing initial q and qDot for ", macroscopic%N, " particles, with tEquiv ", tEquiv
   call init(real(1.2**tEquiv))
   !call init(spooky=1)

   write (*,*) "Done!\n"

   !histOutput = hist(qDot(2,1:N),10)
   !write(*,*) "THIS IS WHAT I SENT: ",qDot(2,1:N)
   !call nextPlot2d(histOutput(1:10,1),histOutput(1:10,2))

   !write(*,*) qDot(2,1:N)

   call setXrange(0.0,real(boxSize))
   call setYrange(0.0,real(boxSize))
   call setZrange(0.0,real(boxSize))

   write (*,*) "Starting iterations [will do ", tN," of em] (simulating ",dt*tN," physical seconds)..."
   !lastFrame=frame
   call cpu_time(temp1)   
   !collisionAlgorithm
   call iterateTheseManyTimes(tN,plotGraphs=1)
   call cpu_time(temp2)
   write(*,*) "Done iterating in  ", temp2- temp1, " seconds\n"
   !write (*,*) "Done iterating :) \n"
   !l=l+1
   
      

   !SHUT THIS FOR SOMETIME
   ! macroscopic%P=sum(framesRecord(10:tN)%P)/real(tN-10)
   ! !volume is known
   ! macroscopic%E=sum(framesRecord(10:tN)%E)/real(tN-10)   
   ! !avgE=sum(E(1000:tN))/real(tN-1000)
   ! !temp1=avgP*volume - (N*Kb*temperature(avgE))
   
   ! temp1=macroscopic%P*macroscopic%V - (2/3.0)*macroscopic%E
   ! write(4,*) macroscopic%N,temp1,macroscopic%P,macroscopic%E,macroscopic%P*macroscopic%V,(2.0/3.0)*macroscopic%E




   !windowSize(l)=evaluateAvgWindowSizeForP(Pw,15.0)
   !windowSizeN(l)=N

   !histOutput = hist(qDot(2,1:N),10)
   !call nextPlot2d(histOutput(1:10,1),histOutput(1:10,2))

   !write(*,*) q(2,1:N)
   !E(1:tN)
   
   !histOutput=hist(framesRecord(1:tN)%E,50)
   
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

! call plot2dSave(framesRecord(1:tN)%t,framesRecord(1:tN)%qT(3),"qT3.pdf",picFormat=1)
! call plot2dSave(framesRecord(1:tN)%t,framesRecord(1:tN)%qT(2),"qT2.pdf",picFormat=1)
! call plot2dSave(framesRecord(1:tN)%t,framesRecord(1:tN)%qT(1),"qT1.pdf",picFormat=1)


!  call plot2dSave(t,qT(2,:),"qT2.pdf")
!  call plot2dSave(t,qT(3,:),"qT3.pdf")
!  call plot2dSave(t,E,"energyT.pdf")
! call plot2dSave(framesRecord(1:tN)%t,framesRecord(1:tN)%qDotT(3),"qDotT3.pdf",picFormat=1)
! call plot2dSave(framesRecord(1:tN)%t,framesRecord(1:tN)%qDotT(2),"qDotT2.pdf",picFormat=1)
! call plot2dSave(framesRecord(1:tN)%t,framesRecord(1:tN)%qDotT(1),"qDotT1.pdf",picFormat=1)



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
call system("xdg-open result3d.avi")

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
  subroutine init(parametricT,spooky)
    real, optional:: parametricT
    integer, optional :: spooky
    !intialize boxSize etc. required for iterating
    boxSize=macroscopic%V**(1.0/3.0)
    area=macroscopic%V**(2.0/3.0)
    
    !box size is between 0 and 1
    !call srand(1000)
    
    !initialize particles to be at some random location
    framesRecord(1)%E=0
    framesRecord(:)%t=0

    do i=1,macroscopic%N
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


       if(present(spooky)) then
       !Its the spooknfiguration
          frame%particles(i)%q(1)=((real(i)/real(macroscopic%N))*(boxSize-(2*radius))) + radius
          frame%particles(i)%q(2)=2*radius
          frame%particles(i)%q(3)=2*radius

          frame%particles(i)%qDot(:)=0.0
          frame%particles(i)%qDot(1)=-5.0 !abs(randomNormal()*10.0)
          
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
    integer(kind=4) :: ii,s
    !this is for efficient collision detection
    integer(kind=4), dimension(3) :: qDisc
    real, dimension(3) :: qVector,qDotVector,unitqVector,qDotNewI,qDotNewII
    logical :: collision
    !the a,b,c are for solving a quadratic! heheh
    real :: deltaT,a,b,c
    !this is for neighbour searching
    !integer, parameter :: nPrec=20 !,density=10
    integer :: nPrec
    integer :: allocateStatus
    integer, dimension(:,:,:), allocatable :: nMat
    !integer(kind=4), dimension(0:nPrec,0:nPrec,0:nPrec,0:density)::nMat=0
    integer  :: nMatCurr=0
    integer :: rk
    integer(kind=4) :: densityInBox,ni
    type (particleLike) :: tempParticle
    type temporary
       real, dimension(3) :: k1,k2,k3,k4
    end type temporary
    integer :: lastFrameNumber = 0
    type(temporary) :: temp
    integer :: neighbourCount
    integer, dimension(100) :: neighbourAddress
    integer, dimension(:,:), allocatable::boxAddresses
    integer, dimension(3):: cAdd
    !call nMatPrettyPrint(nMat)
    call initProgress()
    !write(*,'(A)') "\b[##-------]"
    !to have roughly as many grid cells as the number of particles
    nPrec=(macroscopic%N)**(1.0/3)
    allocate(nMat(0:nPrec,0:nPrec,0:nPrec),stat=allocateStatus)    
    if(allocateStatus /= 0) stop "Not enough memory :("

    time=0
    !lastFrame=frame
    do k=1,numberOfIterations
       framesRecord(k)%t=time
       time=time+dt
       framesRecord(k)%E=0
       framesRecord(k)%P=0
       frame%pressureWall(:,:)=0
       lastFrame%pressureWall(:,:)=0
       !collision with the walls
       !initialize neigbhour matrix to zero
       nMat=0
       !these two innoscent looking buggers are the bottleneck for the speed!
       !frame%particles(:)%next=0
       !frame%particles(:)%collided=0
       !write(*,*) nMat

       !populates the grid for nearest neighbour search
       do i=1,macroscopic%N
          lastFrame%particles(i)=frame%particles(i)

          !initialize for each particle some quantities
          frame%particles(i)%next=0
          frame%particles(i)%collided=0

          !First discritize the particle's location
          qDisc=int((lastFrame%particles(i)%q(:)/boxSize)*nPrec)
          !extract the relavent box from the neighbour matrix
          if( inRange(qDisc,0,nPrec)) then
             nMatCurr=nMat(qDisc(1),qDisc(2),qDisc(3))
          else
             stop "Particle went out of the box, consider decreasing dt"
          end if

          ! if(nMatCurr > macroscopic%N) then
          !    write(*,*) nMat(qDisc(1),qDisc(2),qDisc(3))
          !    stop "Initialization issue?"
          ! end if
          
          !IF the box was empty, coolio :)
          if(nMatCurr==0) then
             !store the address of the first particle!
             nMatCurr=i
             ! if(nMatCurr>macroscopic%N) then
             !    stop "why is i greater than N"
             ! end if
          else
             !we loop over all particles inside the box
             !since nMatCurr is not zero, (you checked it), so the number in nMatCurr must be the address/index of the first particle in the box (by construction)                   
             !let ii hold that address/index (of the first particle in the box)                   
             ii=nMatCurr                   
             do while (ii .ne. 0)
                tempParticle=frame%particles(ii)
                !Find the next particle                      
                if(tempParticle%next .ne. 0) then                         
                   !update the index to point to the next particle
                   ii=tempParticle%next
                else
                   !ii holds the index to the last particle
                   !to it i'll add my ith particle (one i'm looking at)
                   frame%particles(ii)%next=i
                   !ii=0 !(because there's no next particle!)
                   exit                         
                end if
             end do
          end if
          nMat(qDisc(1),qDisc(2),qDisc(3))=nMatCurr

       end do
       
       do rk=1,4
          do i=1,macroscopic%N
             !save the old positions in the last frame
             !lastFrame%particles(i) = frame%particles(i)

             !GET A LIST OF NEIGHBOURS          
             neighbourCount=0

             !First discritize the particle's location
             qDisc=int((lastFrame%particles(i)%q(:)/boxSize)*nPrec)
             !get a list of neibhour boxes depending on range of interaction
             boxAddresses=getBoxAddresses(4*radius,boxSize/real(nPrec))
             !write(*,*) size(boxAddresses(:,1))
             do j=1,size(boxAddresses(:,1))
                !goto the neibhouring box
                cAdd=qDisc + boxAddresses(j,:)
                !if it exists ofcourse
                if(cAdd(1)>=0 .and. cAdd(2)>=0 .and. cAdd(3)>=0 .and. cAdd(1)<=nPrec .and. cAdd(2)<=nPrec .and. cAdd(3)<=nPrec) then
                   !extract the relavent box from the neighbour matrix          
                   nMatCurr=nMat(cAdd(1),cAdd(2),cAdd(3))
                   !write(*,*) nMatCurr
                   if(nMatCurr .ne. 0) then
                      neighbourCount=neighbourCount+1
                      neighbourAddress(neighbourCount)=nMatCurr
                      ii=nMatCurr
                      do while (ii .ne. 0)
                         tempParticle=frame%particles(ii)
                         if(tempParticle%next .ne. 0) then                         
                            neighbourCount=neighbourCount+1
                            neighbourAddress(neighbourCount)=tempParticle%next
                            ii=tempParticle%next
                         else
                            exit                         
                         end if
                      end do
                   end if
                   !neighbourAddress
                end if
             end do

             ! write(*,*) " for ",i, neighbourAddress(1:neighbourCount)
             ! write(*,*) " for ",i, qDisc
             ! write(*,*) " for ",i, size(boxAddresses(:,1))

             !firs set the force being applied on the ith particle to be zero
             frame%particles(i)%force=0
             !update particle velocity (bad algorithm)
             !to do that, first find the force being applied on the particle

             !EULER
             ! frame%particles(i)%force = effectiveForce(lastFrame%particles(i)%q,i,lastFrame%particles, (/ (s, s=1,macroscopic%N) /), potentialParameter)

             !RK 4
             if(rk==1) then
                frame%particles(i)%k1=frame%particles(i)%qDot
                !the second parameter in effectiveF gives which q to use for evaluating the force
                frame%particles(i)%m1=effectiveF(i,1,frame%particles,neighbourAddress(1:neighbourCount),potentialParameter)

                frame%particles(i)%q2=frame%particles(i)%q + 0.5*dt*frame%particles(i)%k1
                !the following in general is needed, but in this case, it is not (since the force doesn't depend on the velocity)
                !frame%particles(i)%qDot2=frame%particles(i)%qDot + 0.5*dt*frame%particles(i)%m1
             else if(rk==2) then
                frame%particles(i)%k2=frame%particles(i)%qDot + 0.5*dt*frame%particles(i)%m1
                frame%particles(i)%m2=effectiveF(i,2,frame%particles,neighbourAddress(1:neighbourCount),potentialParameter)

                frame%particles(i)%q3=frame%particles(i)%q + 0.5*dt*frame%particles(i)%k2
             else if(rk==3) then
                frame%particles(i)%k3=frame%particles(i)%qDot + 0.5*dt*frame%particles(i)%m2
                frame%particles(i)%m3=effectiveF(i,3,frame%particles,neighbourAddress(1:neighbourCount),potentialParameter)

                frame%particles(i)%q4=frame%particles(i)%q + dt*frame%particles(i)%k3
             else if(rk==4) then
                frame%particles(i)%k4=frame%particles(i)%qDot + dt*frame%particles(i)%m3
                frame%particles(i)%m4=effectiveF(i,4,frame%particles,neighbourAddress(1:neighbourCount),potentialParameter)

                frame%particles(i)%q=frame%particles(i)%q + (frame%particles(i)%k1 + (2*frame%particles(i)%k2) + (2*frame%particles(i)%k3) + frame%particles(i)%k4)*dt/6.0
                frame%particles(i)%qDot=frame%particles(i)%qDot + (frame%particles(i)%m1 + (2*frame%particles(i)%m2) + (2*frame%particles(i)%m3) + frame%particles(i)%m4)*dt/6.0
                !frame%particles(i)%q4=frame%particles(i)%q + dt*frame%particles(i)%k3


                !AND NOW do what you want to about the wall collisions etc.
                !check for wall collisions and evaluate pressure on each wall on the fly
                do j=1,3
                   ! if (q(j,i) >= 1.0) then
                   ! if (q(j,i) >= boxSize) then
                   if (frame%particles(i)%q(j)>=boxSize-radius) then
                      frame%pressureWall(j,1)=frame%pressureWall(j,1)+ ( (2*m*frame%particles(i)%qDot(j))/ (dt*area) )
                      frame%particles(i)%qDot(j)=-frame%particles(i)%qDot(j)
                      frame%particles(i)%q(j)=frame%particles(i)%q(j)-2*(frame%particles(i)%q(j)-(boxSize-radius))
                      if(frame%particles(i)%q(j)<=radius) then
                         !stop "Decrease dt! After reflecting from the wall and correcting, the particle is going outside the box! :("
                      end if
                      !else if (q(j,i) <= 0.0) then
                   else if (frame%particles(i)%q(j)<=radius) then
                      frame%pressureWall(j,2)=frame%pressureWall(j,2)-( (2*m*frame%particles(i)%qDot(j)) / (dt*area) )
                      frame%particles(i)%qDot(j)=-frame%particles(i)%qDot(j)
                      frame%particles(i)%q(j) = (2*radius)-frame%particles(i)%q(j) 
                      if(frame%particles(i)%q(j)>=boxSize-radius) then
                         !stop "Decrease dt! After reflecting from the wall and correcting, the particle is going outside the box! :("
                      end if
                   end if
                end do
                !evaluate total energy
                framesRecord(k)%E=framesRecord(k)%E+energy(frame%particles(i)%qDot(:))

             end if
                
                
             ! temp%k1 = (effectiveForce(lastFrame%particles(i)%q,i,lastFrame%particles, neighbourAddress, potentialParameter, neighbourCount)/m)*dt
             ! temp%k2 = (effectiveForce(lastFrame%particles(i)%q + 0.5*temp%k1,i,lastFrame%particles, neighbourAddress, potentialParameter, neighbourCount)/m)*dt
             ! temp%k3 = (effectiveForce(lastFrame%particles(i)%q + 0.5*temp%k2,i,lastFrame%particles, neighbourAddress, potentialParameter, neighbourCount)/m)*dt
             ! temp%k4 = (effectiveForce(lastFrame%particles(i)%q + temp%k3,i,lastFrame%particles, neighbourAddress, potentialParameter, neighbourCount)/m)*dt

             ! temp%k1 = (effectiveForce(lastFrame%particles(i)%q,i,lastFrame%particles, (/ (s, s=1,macroscopic%N) /), potentialParameter)/m)*dt
             ! temp%k2 = (effectiveForce(lastFrame%particles(i)%q + 0.5*temp%k1,i,lastFrame%particles, (/ (s, s=1,macroscopic%N) /), potentialParameter)/m)*dt
             ! temp%k3 = (effectiveForce(lastFrame%particles(i)%q + 0.5*temp%k2,i,lastFrame%particles, (/ (s, s=1,macroscopic%N) /), potentialParameter)/m)*dt
             ! temp%k4 = (effectiveForce(lastFrame%particles(i)%q + temp%k3,i,lastFrame%particles, (/ (s, s=1,macroscopic%N) /), potentialParameter)/m)*dt

             !OBSOLETE
             ! do s=1,macroscopic%N
             !    if(s .ne. i) then
             !       frame%particles(i)%force=frame%particles(i)%force + force(lastFrame%particles(s)%q,lastFrame%particles(i)%q,potentialParameter)
             !    end if
             ! end do

             ! now with this, find the particle velocity
             !EULER
             ! frame%particles(i)%qDot=frame%particles(i)%qDot + ((frame%particles(i)%force/m)*dt)


             !RK 4
             !frame%particles(i)%qDot=frame%particles(i)%qDot + (1/6.0)*(temp%k1 + 2*temp%k2 + 2*temp%k3 + temp%k4)

             ! update particle position
             ! part of RK 4, but euler
             ! frame%particles(i)%q(:)=frame%particles(i)%q(:)+(lastFrame%particles(i)%qDot(:)*dt)

             !EULER (MIXED)
             !frame%particles(i)%q(:)=frame%particles(i)%q(:)+(frame%particles(i)%qDot(:)*dt)



          end do
       end do
       

       !THis is to save the location and position of particle 1
       !at different times
       framesRecord(k)%qT(:)=frame%particles(1)%q(:)
       framesRecord(k)%qDotT(:)=frame%particles(1)%qDot(:)
       !average pressure over all walls (or equivalently, could've put the surface area as 6*area of square)
       framesRecord(k)%P=sum(sum(frame%pressureWall,dim=1),dim=1)/6.0
       if (present(plotGraphs)) then
         if(time<renderDuration) then 
             if(plotGraphs==1 .or. plotGraphs==3) then
                if ( int(time*fps) .ne. lastFrameNumber ) then
                   call nextPlot3d(frame%particles(1:macroscopic%N)%q(1),frame%particles(1:macroscopic%N)%q(2),frame%particles(1:macroscopic%N)%q(3))!q(2,1:N),q(3,1:N))
                   lastFrameNumber=int(time*fps)
                end if
             end if
          end if
          if(k<mN) then
             if(plotGraphs==2 .or. plotGraphs==3) then
                
                histOutput = hist((frame%particles(1:macroscopic%N)%qDot(1)**2) + (frame%particles(1:macroscopic%N)%qDot(2)**2) + (frame%particles(1:macroscopic%N)%qDot(3)**2),20)
                ! hist((qDot(2,1:N)*qDot(2,1:N))+ (qDot(1,1:N)*qDot(1,1:N)) + (qDot(3,1:N)*qDot(3,1:N)),20)
                call nextPlot2d(histOutput(1:20,1),histOutput(1:20,2))
             end if
          end if
       end if
       !write(*,'(f5.5A)',advance="no") real(k)/real(numberOfIterations),"\r"
       call incProgress(k,numberOfIterations)
       !lastFrame=frame
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

  function V(r,q,par)
    type(potentialParametersLike) ::par
    real, dimension(3):: r,q
    real :: modr,V
    modr=sqrt(sum((r-q)**2))
    !if( sum(r*r - q*q) < d*d) then
    V=par%epsilon * ( (par%rm/modr)**12 - 2*(par%rm/modr)**6)       
  end function V

  function Vprime(r,q,par)
    type(potentialParametersLike) ::par
    real, dimension(3):: r,q
    real :: modr,Vprime
    modr=sqrt(sum((r-q)**2))
    ! Vprime=-12.0* (par%epsilon) *(1/modr) * ( ((par%rm/modr)**12) - ((par%rm/modr)**6))
    Vprime=-(par%epsilon) *(1/modr) * ( par%X*par%A*((par%rm/modr)**par%A) - par%Y*par%B*((par%rm/modr)**par%B))
  end function Vprime


  
  function force(r,q,par)
    type(potentialParametersLike) ::par
    real, dimension(3):: r,q,force,rCap
    real :: modr
    rCap=(q-r) /sqrt( sum( (q-r)**2 ))
    force=-Vprime(r,q,par) * rCap
  end function force


  function effectiveForce(q,i,particles,particleIDs,par,numParticles)
    real, dimension(3),intent(in) :: q
    real, dimension(3):: effectiveForce,r
    integer, intent(in):: i
    integer :: s,num
    integer, optional :: numParticles
    type(particleLike), dimension(:),intent(in) :: particles
    integer, dimension(:),intent(in) :: particleIDs
    type(potentialParametersLike),intent(in) :: par
    num=size(particleIDs)
    if(present(numParticles)) then
       num=numParticles
    end if
    effectiveForce=0
    do s=1, num       
       if(particleIDs(s) .ne. i) then
          r=particles(particleIDs(s))%q
          !effectiveForce=effectiveForce + force(particles(particleIDs(s))%q,q,par)
          effectiveForce=effectiveForce + force(r,q,par)
          !write(*,*) particles(particleIDs(s))%q
          if(particleIDs(s) > macroscopic%N) then
             write(*,*) particleIDs(s)
          end if
          !effectiveForce=effectiveForce + (/ 0.0,0.0,0.0 /)
       end if
    end do
  end function effectiveForce

  function effectiveF(i,j,particles,particleIDs,par,numParticles)
    real, dimension(3):: q
    real, dimension(3):: effectiveF,r
    integer, intent(in):: i,j
    integer :: s,num
    integer, optional :: numParticles
    type(particleLike), dimension(:),intent(in) :: particles
    integer, dimension(:),intent(in) :: particleIDs
    type(potentialParametersLike),intent(in) :: par
    num=size(particleIDs)
    if(present(numParticles)) then
       num=numParticles
    end if
    effectiveF=0
    do s=1, num       
       if(particleIDs(s) .ne. i) then
          if(j==1) then
             r=particles(particleIDs(s))%q
             q=particles(i)%q
          else if(j==2) then
             r=particles(particleIDs(s))%q2
             q=particles(i)%q2
          else if(j==3) then
             r=particles(particleIDs(s))%q3
             q=particles(i)%q3
          else if(j==4) then
             r=particles(particleIDs(s))%q4
             q=particles(i)%q4
          end if
          !effectiveForce=effectiveForce + force(particles(particleIDs(s))%q,q,par)
          effectiveF=effectiveF + force(r,q,par)
          !write(*,*) particles(particleIDs(s))%q
          if(particleIDs(s) > macroscopic%N) then
             write(*,*) particleIDs(s)
          end if
          !effectiveForce=effectiveForce + (/ 0.0,0.0,0.0 /)
       end if
    end do
  end function effectiveF


  ! !gradV
  ! function gradV (r,q,par)
    
  ! end function gradV

  function numberOfBoxes(R,B,rrr)
    real:: R,B,s
    integer :: rr,bb,numberOfBoxes
    integer, optional :: rrr
    s=B/2.0
    bb=B/s

    rr=ceiling(R/s)
    !write(*,*) "rr=", R/s,rr
    if(mod(rr,2).eq. 0) then
       rr=rr+1
    end if
    numberOfBoxes=rr**3
    if(present(rrr)) then
       rrr=(rr-1)/2
    end if
  end function numberOfBoxes

  function getBoxAddresses(R,B)
    real:: R,B
    integer:: maxNeBox,maxBoxIndex,i,j,k,m
    integer, dimension(:,:), allocatable :: getBoxAddresses,boxAddresses
    maxNeBox=numberOfBoxes(R,B,maxBoxIndex) !((2*rangeOfForce)/smallBoxSize)**3.0


    ! if( allocated(boxAddresses)) then
    !    deallocate(boxAddresses)
    ! end if
    allocate(boxAddresses(maxNeBox,3))

    ! write(*,*) "For radius ", rangeOfForce, " and smallBoxSize ", smallBoxSize, " I'll need atmost ", maxNeBox, " neighbouring boxes"
    ! write(*,*) "I'll loop from 0 to ", maxBoxIndex

    m=1
    do i=-maxBoxIndex,maxBoxIndex
       do j=-maxBoxIndex,maxBoxIndex
          do k=-maxBoxIndex,maxBoxIndex
             if( (i**2 + j**2 + k**2) <= (maxBoxIndex)**2 ) then
                !write(*,*) i,j,k
                boxAddresses(m,1)=i
                boxAddresses(m,2)=j
                boxAddresses(m,3)=k
                m=m+1
             end if
          end do
       end do
    end do

    if(allocated(getBoxAddresses)) then
       deallocate(getBoxAddresses)
    end if
    allocate(getBoxAddresses(m-1,3))
    getBoxAddresses(:,:)=boxAddresses(1:m-1,:)
    
  end function getBoxAddresses


  function inRange(vector,min,max)
    integer(kind=4), dimension(3):: vector
    integer:: min, max
    logical :: inRange
    inRange=.false.
    if (min<=vector(1) .and. vector(1)<=max .and. min<=vector(2) .and. vector(2)<=max .and. min<=vector(3) .and. vector(3)<=max) then
       inRange=.true.
    end if
  end function inRange
  
end program bouncer


!make x array
!dx = (x_end - x_start)/n        
!x(0:n) = [(i*dx, i=0,n)]        
