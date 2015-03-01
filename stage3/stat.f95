module stat
  implicit none
  
contains  
  !hist(:,1) this is the value
  !hist(:,2) this is the frequency
  function hist(data,numberOfBins)
    real, dimension(:), intent(in) :: data
    integer, intent(in) :: numberOfBins
    integer :: i
    real :: binSize, left, right, dataMax, dataMin
    real, dimension(numberOfBins,2) :: hist
    !real, dimension(10000,2) :: hist
    dataMax = maxval(data)
    dataMin = minval(data)
    !write(*,*) "THIS IS WHAT GOT SENT",data
    !write(*,*) dataMax,dataMin
    if(dataMax==dataMin) then
       dataMax=dataMax+0.1
       dataMin=dataMin-0.1
    end if
    binSize = (dataMax - dataMin)/real(numberOfBins)    
    ! hist(1,1)=0
    ! hist(1,2)=0
    !put i=1
    do i = 1, numberOfBins
       left=dataMin+((i-1)*binSize)
       right=dataMin + (i*binSize)
       hist(i,2) = count( left<=data .and. data<right )
       hist(i,1) = (left + right)/2
    end do
  end function hist

  function stdDevFromHist(hist)
    real :: expX,expXsquare,stdDevFromHist
    real, intent(in), dimension(:,:) :: hist
    !find expectation
    expX=sum(hist(:,1)*hist(:,2))/sum(hist(:,2))
    expXsquare=sum(hist(:,1)*hist(:,1)*hist(:,2))/sum(hist(:,2))
    stdDevFromHist=sqrt(expXsquare - expX*expX)    
  end function stdDevFromHist

  function randomNormal()
    real :: tempSum, randomNormal
    integer :: i,iMax=1000
    tempSum=0
    do i=1,iMax
       tempSum=tempSum + rand()
    end do
    !to get it centred at the origin
    randomNormal=(tempSum/real(iMax)-0.5)
  end function randomNormal
  
end module stat
