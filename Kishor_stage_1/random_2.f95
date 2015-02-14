program rand_2
	implicit none
	integer:: die
	real:: x
	integer i,j(1)
	call system_clock(i)
	j(1)=i
	call random_seed(put=j)
	call random_number(x)
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
end program rand_2
