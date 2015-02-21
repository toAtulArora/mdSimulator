module stat
  implicit none
  
contains  
  function hist(data,numberOfBins)
    real, dimension(:), intent(in) :: data
    integer, intent(in) :: numberOfBins
    integer :: dataMax, dataMin, i
    real :: binSize, left, right
    real, dimension(numberOfBins,2) :: hist
    !real, dimension(10000,2) :: hist
    dataMax = maxval(data)
    dataMin = minval(data)
    binSize = (dataMax - dataMin)/real(numberOfBins)    
    do i = 1, numberOfBins
       left=dataMin+((i-1)*binSize)
       right=dataMin + (i*binSize)
       hist(i,2) = count( left<=data .and. data<right )
       hist(i,1) = (left + right)/2
    end do
  end function hist
  
end module stat
