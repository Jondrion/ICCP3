module PLplot3d
  
use plplot
implicit none

contains

subroutine plot_init()
    real(kind=plflt) arrow2_x(6),arrow2_y(6)

    data arrow2_x/-0.5_plflt, 0.3_plflt, 0.3_plflt, 0.5_plflt, 0.3_plflt, 0.3_plflt/
    data arrow2_y/0._plflt, 0._plflt, 0.2_plflt, 0._plflt, -0.2_plflt, 0._plflt/


    call plsdev("xcairo")
    call plinit()
    
    call plsvect(arrow2_x, arrow2_y, .false.)

    call plspause(.false.)
    !call plend()

end subroutine plot_init

subroutine plot_points(vector, x, y)
    integer, intent(in) :: x,y
    real(8), intent(in) :: vector(y,x,2)
    real(8) :: xx, yy
    real(kind=plflt) xmin, xmax, ymin, ymax;
    real(kind=plflt) xg(x,y),yg(x,y),u(x,y),v(x,y)
    integer :: i,j,k
    real(kind=plflt) :: clev(11)
    character(len=1)  defined



    do i=1,x
        xx = i
        do j=1,y
          yy = j
          xg(i,j) = xx
          yg(i,j) = yy
          !u(i,j) = yy
          !v(i,j) = -xx
        enddo
    enddo

    xmin=0
    xmax=x+1
    ymin=0
    ymax=y+1

    do k=1,11
        clev(k)=minval(u)+(k-1)/10*(maxval(u)-minval(u))
    enddo 

    call plclear()
    call plcol0(2)
    call plenv(xmin, xmax, ymin, ymax, 0, 0)    
    call plcol0(1)
    call plvect(vector(:,:,1),vector(:,:,2),0.0_plflt,xg,yg)
    !call plshades(u, defined, xmin, xmax, ymin, ymax, clev, 0._plflt, 1, 1._plflt, .false.)
        
    
    
    call plflush()
end subroutine

end module PLplot3d
