module PLplot3d
  
use plplot
implicit none

contains

subroutine plot_init()
    integer, parameter :: nc=15
    integer :: col1(nc+1), col2(nc+1), col3(nc+1),k
    real(kind=plflt) arrow2_x(6),arrow2_y(6)

    data arrow2_x/-0.5_plflt, 0.3_plflt, 0.3_plflt, 0.5_plflt, 0.3_plflt, 0.3_plflt/
    data arrow2_y/0._plflt, 0._plflt, 0.2_plflt, 0._plflt, -0.2_plflt, 0._plflt/

    call plscolbg(255,255,255)
    call plscol0(7,0,0,0)


    call plsdev("xcairo")
    call plinit()
    
    call plsvect(arrow2_x, arrow2_y, .false.)

    ! -- change color sheme for plshades to blue and orange
    do k=1,nc+1
       col1(k)=int(dble(nc+1-k)/dble(nc)*255)
       col2(k)=100
       col3(k)=int(dble(k-1)/dble(nc)*255)
       call plscmap1( col1 ,col2 ,col3)
    enddo

    call plspause(.false.)
    !call plend()

end subroutine plot_init

subroutine plot_points(vector, x, y, X_object, n_vertices, CoM)

use dispmodule
    integer, intent(in) :: x,y, n_vertices
    real(kind=plflt), intent(in) :: vector(y,x,2), X_object(n_vertices,2), CoM(2)
    real(8) :: xx, yy
    real(kind=plflt) xmin, xmax, ymin, ymax;
    real(kind=plflt) xg(y,x),yg(y,x),line(x)
    integer :: i,j,k
    integer, parameter :: nc=15
    real(kind=plflt) clev(nc+1)

    character(len=1)  defined


    do i=1,x
        xx = i
        line(i)=y*sqrt(3._8)/2._8
        do j=1,y
          yy = j
          xg(j,i) = xx+modulo(j+1,2)*0.5
          yg(j,i) = yy*sqrt(3._8)/2._8
        end do
    end do

    xmin=0._plflt
    xmax=xg(2,x)+1._plflt
    ymin=0._plflt
    ymax=yg(y,1)+0._plflt

    do k=1,nc+1
       clev(k)=minval(vector(2:y-1,:,1))+dble(k-1)/dble(nc)*(maxval(vector(2:y-1,:,1))-minval(vector(2:y-1,:,1)))
    enddo 

    !call plclear()

    call plcol0(7)
    call plenv(xmin, xmax, ymin, ymax, 1, 0)
    call PLshades(vector(2:y-1,:,1), defined, xmin, xmax, ymin, ymax, clev,1._plflt, 1, 0._plflt,xg(2:y-1,:),yg(2:y-1,:)) 

    call plcol0(3)
    call plvect(vector(:,:,1),vector(:,:,2),8.0_plflt,xg,yg)

    call plcol0(7)
    !call plline(xg(1,:),line-1 )
    !call plline(xg(1,:),line-y+2)

    call plfill( X_object(:,1), X_object(:,2) )
    call plcol0(1)
    call plpoin([CoM(1),CoM(1)],[CoM(2),CoM(2)],1)


    
    call plflush()
end subroutine

end module PLplot3d
