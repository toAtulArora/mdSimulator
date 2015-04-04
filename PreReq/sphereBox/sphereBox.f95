program sphereBox
  implicit none
  integer, parameter :: nPrec=100
  real, parameter :: rangeOfForce=3,smallBoxSize=1.0
  integer, dimension(nPrec,nPrec,nPrec):: nMat
  integer, dimension(:,:), allocatable :: boxAddresses
  !maximum neignbouring boxes
  integer :: k
  boxAddresses=getBoxAddresses(rangeOfForce,smallBoxSize)
  ! write(*,*) boxAddresses
  do k=1, size(boxAddresses(:,1))
     write(*,*) boxAddresses(k,:)
  end do
  !write (*,*) boxAddresses(1:m,:)
  !list boxList
contains
  ! subroutine boxList(boxSize,radius)
  !   write(*,*) "I was called"
  ! end subroutine boxList

  !R is the range of interaction
  !B is the small box size (grid size)
  !rrr is the range you'd iterate to: 0 to rrr
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
  
  
  ! function isBoxInRange(R,B,rx)
  !   real:: R,B,s
  !   real, dimension(3)::Rx
  !   integer,dimension(3) :: rx
  !   logical :: isBoxInRange = .false.
  !   s=B/2.0
  !   Rx=s*rx
    
  ! end function isBoxInRange



  
end program sphereBox
