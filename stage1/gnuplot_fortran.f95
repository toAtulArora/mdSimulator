module gnuplot_fortran
  implicit none
  integer :: graphCount=0
  character(13) :: fileName
contains
  subroutine plot2d(x,y,filename)
    real, intent(in), dimension(:) :: x,y
    character, intent(in), dimension(13) :: filename
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
    write(2,*) "set terminal jpeg"
    write(2,*) "set output 'temp/",filename,"'"
    write(*,*) "The filename you gave was: ", filename
    write(2,*) "plot 'tempData.dat'"
    close(2)

    
    call system ("gnuplot 'command'")
    call system ("rm tempData.dat")
  end subroutine plot2d

  subroutine plot3d(x,y,z,filename)
    real, intent(in), dimension(:) :: x,y,z
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
    write(*,*) "The filename you gave was: ", filename
    write(2,*) "set xrange [-1:10]"
    write(2,*) "set yrange [-1:10]"
    write(2,*) "set zrange [-1:10]"
    write(2,*) "splot 'tempData.dat' using 1:2:3"
    close(2)

    
    call system ("gnuplot 'command'")
    call system ("rm tempData.dat")
  end subroutine plot3d

  subroutine nextPlot2d(x,y)
    real, intent(in), dimension(:) :: x,y
    graphCount = graphCount+1
    write(fileName,'(a,i4.4,a)') 'file',graphCount,'.jpeg'
    call plot2d(x,y,fileName)
  end subroutine nextPlot2d

  subroutine nextPlot3d(x,y,z)
    real, intent(in), dimension(:) :: x,y,z
    graphCount = graphCount+1
    write(fileName,'(a,i4.4,a)') 'file',graphCount,'.jpeg'
    call plot3d(x,y,z,fileName)
  end subroutine nextPlot3d

  subroutine startPlot()
    call system ("rm -r temp")
    call system ("mkdir temp")    
    graphCount=0
  end subroutine startPlot
  subroutine endPlot()
    write (*,*) "Converting to avi.."
    call system ("avconv -i 'temp/file%04d.jpeg' result.avi")
    write (*,*) "Done!"
    !TODO: command for creating a movie and cleanup
  end subroutine endPlot
!  subroutine createPlot()

end module gnuplot_fortran
