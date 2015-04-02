program typeD
  type coordinate
     real, dimension(3):: x
  end type coordinate

  type particle
     type(coordinate) q,qDot
     type(particle), pointer :: next
  end type particle

  ! real, target :: abc=100
  ! real, pointer :: pAbc

  type(coordinate) qDef!( (/ 0.0,1.0,2.0 /) )
  type(particle) Chanchuon
  type(particle), dimension(0:10), target :: Chanchuons
  type(particle), pointer :: temp
  integer :: k



  ! write(*,*)" input coordinate x"
  ! read(*,*) Chanchuon%q
  
  ! write(*,*)"input velocity: z"
  ! read(*,*) Chanchuon%qDot

  ! write(*,*)"printing the object", Chanchuon

  ! write(*,*) "i can access q " ,Chanchuon%q
  !allocate (abc)
  ! pAbc => abc
  ! abc=100
  ! write(*,*) pAbc
  
  qDef%x= (/ 5.0,1.0,2.0 /)
  do k=0,10
     Chanchuons(k)%q=qDef
     Chanchuons(k)%qDot%x = (qDef%x+qDef%x)
     write(*,*) Chanchuons(k)%q, Chanchuons(k)%qDot
     if(k<10) then
        allocate(Chanchuons(k)%next )
        Chanchuons(k)%next  => Chanchuons(k+1)
     end if
     !write(*,*) "one" , ChanchuonNew(Chanchuon)
     !write(*,*) "two", Chanchuon
  end do
  Chanchuons(k)%next => null()
  Chanchuons=ChanchuonNew(Chanchuons)
  do k=0,10
     write(*,*) Chanchuons(k)%q, Chanchuons(k)%qDot
     ! write(*,*) Chanchuons(k)   
  end do


  ! write(*,*) (Chanchuons(k),k=1,10)
contains

  elemental function ChanchuonNew(Chanchuon)
    type(particle), target, intent(in) :: Chanchuon
    type(particle), target :: ChanchuonNew
    type(particle), target :: tempChancha
    !ChanchuonNew=Chanchuon
    !tempChancha = 
    ChanchuonNew%next => null() !Chanchuon%next
    ChanchuonNew%q=Chanchuon%qDot
    ChanchuonNew%qDot=Chanchuon%q
    !write(*,*) "printing the object"

  end function ChanchuonNew



end program typeD
