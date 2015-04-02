module obj_class
  implicit none
  private
  public :: obj
  integer :: i,j,l=5
  real :: s = 25

  type obj
     integer :: m
     real, dimension (:,:), allocatable :: a
   contains
     procedure :: change_obj
     procedure :: del_obj
  end type obj

  interface obj
     module procedure init_obj
  end interface obj

contains
  type(obj) function init_obj(n)
    integer, intent(in) :: n
    allocate (init_obj%a(n,n))
    init_obj%a=s
    init_obj%m=l
  end function init_obj

  subroutine change_obj(self, n, r)
    class (obj), intent(inout) :: self
    integer, intent(in) :: n
    real, intent(in) :: r
    self%a = r
    self%m = n
  end subroutine change_obj

  subroutine del_obj(self)
    class (obj) :: self
    deallocate(self%a)    
  end subroutine del_obj
end module obj_class
