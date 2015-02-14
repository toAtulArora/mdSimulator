program rand_1	
	implicit none
	integer:: die
	real:: x
	call random_number(x) ! gives a random number in [0,1)
	WRITE (*,*)'RANDOM NUMBER :', x
	die=6*x
	select case (die)
		case (1:3)
			write(*,*) 'Stay in bed'
		case (4,5)
			write(*,*)'Skip lecture'
		case (6)
			write(*,*)'Accuse linguist of laziness'
		case default ! This should never happen
			write(*,*)'Attend lecture'
	end select
end program rand_1
