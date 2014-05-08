subroutine make_object(n_vertices, CoM, object_size, offset, X_object)
    integer, intent(in) :: n_vertices
    real(8), intent(in) :: offset(2), object_size
    real(8), intent(out) :: X_object(n_vertices,2)
    real(8), intent(inout) :: CoM(2)
    real(8) :: Pi
    integer :: i

    Pi=4._8*atan(1d0)



    do i=1, n_vertices
        X_object(i,1)=CoM(1)+object_size*cos(2*Pi*dble(i)/dble(n_vertices))
        X_object(i,2)=CoM(2)+object_size*sin(2*Pi*dble(i)/dble(n_vertices))
    end do

    CoM=CoM+offset

end subroutine