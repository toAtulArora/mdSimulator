module gnuplot_fortran
  implicit none
  integer :: graphCount=0
  character(13) :: fileName
  real, dimension (3) :: rangeStart=(/0.0,0.0,0.0/),rangeEnd=(/1.0,1.0,1.0/)
contains
  subroutine setXrange(rangeStartIn,rangeEndIn)
    real, intent(in) :: rangeStartIn,rangeEndIn
    rangeStart(1)=rangeStartIn
    rangeEnd(1)=rangeEndIn
  end subroutine setXrange

  subroutine setYrange(rangeStartIn,rangeEndIn)
    real, intent(in) :: rangeStartIn,rangeEndIn
    rangeStart(2)=rangeStartIn
    rangeEnd(2)=rangeEndIn
  end subroutine setYrange

  subroutine setZrange(rangeStartIn,rangeEndIn)
    real, intent(in) :: rangeStartIn,rangeEndIn
    rangeStart(3)=rangeStartIn
    rangeEnd(3)=rangeEndIn
  end subroutine setZrange


  subroutine plot2dSave(x,y,filename,rangeXstart,rangeXend,rangeYstart,rangeYend)
    real, intent(in), dimension(:) :: x,y
    character(len=*), intent(in),optional :: filename
    real, intent(in), optional :: rangeXstart, rangeXend, rangeYstart, rangeYend
    integer :: size_x, size_y,i
    size_x = size(x)
    size_y = size(y)
    if (size_x /= size_y ) then
       print *,"Array size mismatch"
    else
       open(unit =1 ,file='tempData.dat')
       do i = 1, size(x)
          write(1,*) x(i), y(i)
       end do
       close(1)
    end if
    
    !open a file to write the commands for gnuplot
    open(unit =2,file='command')
    if(present(filename)) then
       write(2,*) "set terminal jpeg"
       write(2,*) "set output 'temp2d/",filename,"'"
    end if

    if (present(rangeXstart) .and. present(rangeXend) ) then
       write(2,*) "set xrange [",rangeXstart,":",rangeXend,"]"
    end if

    if (present(rangeYstart) .and. present(rangeYend) ) then
       write(2,*) "set yrange [",rangeYstart,":",rangeYend,"]"
    end if

    !write(2,*) "set yrange [0:1]"
    
    write(2,*) "plot 'tempData.dat' w lp"
    close(2)
    
    if (present(filename)) then
       call system ("gnuplot 'command'")
       write(*,*) "'",filename,"' has been saved."
    else
       call system ("gnuplot -persist 'command'")
    end if
    !call system ("gnuplot -persist 'command'")
    !call system ("rm tempData.dat")
  end subroutine plot2dSave

  subroutine plot2d(x,y)
    real, intent(in), dimension(:) :: x,y
    
    integer :: size_x, size_y,i
    size_x = size(x)
    size_y = size(y)
    if (size_x /= size_y ) then
       print *,"Array size mismatch"
    else
       open(unit =1 ,file='tempData.dat')
       do i = 1, size(x)
          write(1,*) x(i), y(i)
       end do
       close(1)
    end if
    
    !open a file to write the commands for gnuplot
    open(unit =2,file='command')
    !write(2,*) "set terminal jpeg"
    !write(2,*) "set output 'temp/",filename,"'"
    !write(2,*) "set yrange [0:1]"
    !write(*,*) "The filename you gave was: ", filename
    write(2,*) "plot 'tempData.dat' w lp"
    close(2)

    !call system ("gnuplot 'command'")
    call system ("gnuplot -persist 'command'")
    call system ("rm tempData.dat")
  end subroutine plot2d

  subroutine plot3d(x,y,z,filename,angle)
    real, intent(in), dimension(:) :: x,y,z
    real, intent(in) :: angle
    character, intent(in), dimension(13) :: filename
    integer :: size_x, size_y, size_z, i
    size_x = size(x)
    size_y = size(y)
    size_z = size(z)
    !TODO: add a condition for z also
    if (size_x /= size_y ) then
       print *,"Array size mismatch"
    else
       open(unit =1 ,file='tempData.dat')
       do i = 1, size(x)
          write(1,*) x(i), y(i), z(i)
       end do
       close(1)
    end if
    
    !open a file to write the commands for gnuplot
    open(unit =2,file='command')
    write(2,*) "set terminal jpeg"
    write(2,*) "set output 'temp/",filename,"'"
    !write(*,"(a)",advance="no") "#"
    !, filename
    write(2,*) "set grid ytics lc rgb '#bbbbbb' lw 1 lt 0"
    write(2,*) "set grid xtics lc rgb '#bbbbbb' lw 1 lt 0"
    write(2,*) "set grid ztics lc rgb '#bbbbbb' lw 1 lt 0"
    write(2,*) "set view 60,",angle
    write(2,*) "set xrange [0:",rangeEnd(1),"]"
    write(2,*) "set yrange [0:",rangeEnd(2),"]"
    write(2,*) "set zrange [0:",rangeEnd(3),"]"
    write(2,*) "splot 'tempData.dat' using 1:2:3"
    close(2)

    
    call system ("gnuplot 'command'")
    call system ("rm tempData.dat")
  end subroutine plot3d

  subroutine nextPlot2d(x,y)
    real, intent(in), dimension(:) :: x,y
    graphCount = graphCount+1
    write(fileName,'(a,i4.4,a)') 'file',graphCount,'.jpeg'
    call plot2dSave(x,y,fileName)
  end subroutine nextPlot2d

  subroutine nextPlot3d(x,y,z)
    real, intent(in), dimension(:) :: x,y,z
    graphCount = graphCount+1
    write(fileName,'(a,i4.4,a)') 'file',graphCount,'.jpeg'
    call plot3d(x,y,z,fileName,mod(real(graphcount),360.0) )
  end subroutine nextPlot3d

  subroutine startPlot()
    call system ("rm -r temp")
    call system ("mkdir temp")    
    call system ("rm -r temp2d")
    call system ("mkdir temp2d")
    graphCount=0
  end subroutine startPlot
  subroutine endPlot()
    write (*,*) "Converting to avi.."
    call system ("avconv -i 'temp/file%04d.jpeg' result.avi")
    !write (*,*) "Done!"
    !TODO: command for creating a movie and cleanup
  end subroutine endPlot
!  subroutine createPlot()

  
end module gnuplot_fortran
