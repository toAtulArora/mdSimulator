program project

	implicit none
	real, parameter :: epsil = 0.1 , bound = 10
	real, parameter :: size_ = ((2*bound/0.1) + 1 ) ! bound determines the number of lattice points in one dimension
	!real :: array_x(size_), array_y(size_), array_z(size_) ! creates an empty array of lattice points
	real :: pos(3) ! placeholder for the position vector
	real :: vel(3) ! placeholder for the velocity vector
	real :: E ! Energy
	integer :: t ! Time
	integer :: i 
	integer :: t_max = 1000

	
	!array_x(1) = 0.0
	!array_y(1) = 0.0
	!array_z(1) = 0.0
	
	!do i = 2,j
		!array_x(i) = array_x(i-1) + epsil
		!array_y(i) = array_y(i-1) + epsil
		!array_z(i) = array_z(i-1) + epsil
	!end do


	! Initial position
	call random_number(pos(1))
	pos(1) = pos(1) * 10 + 0.1*epsil
	call random_number(pos(2))
	pos(2) = pos(2) * 10 + 0.1*epsil
	call random_number(pos(3))
	pos(3) = pos(3) * 10 + 0.1*epsil
	
	! initial velocity
	call random_number(vel(1))
	vel(1) = vel(1) + 0.1*epsil
	call random_number(vel(2))
	vel(2) = vel(2) + 0.1*epsil
	call random_number(vel(3))
	vel(3) = vel(3) + 0.1*epsil
	open(20, file = 'bounce.dat')
	do t = 1,t_max
		pos(1) = pos(1) + vel(1) * 1
		pos(2) = pos(2) + vel(2) * 1
		pos(3) = pos(3) + vel(3) * 1
		write(20,*), pos(1),pos(2), pos(3)
		do i = 1,3
			if ( ((pos(i) + vel(i)) .lt. 0.1 .and. vel(i) .lt. 0) .or. ((10 - pos(i)-vel(i)) .lt. 0 .and. vel(i) .gt. 0.0) ) then
				vel(i) = - vel(i)
			end if
		end do  	
	end do
	close(20)

	
	
end program project
