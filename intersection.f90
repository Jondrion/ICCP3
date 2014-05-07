!-- p1 and p2 are lattice points (p1 is fluid point and p2 is wall point)
!-- p3 and p4 are corner points of wall object
subroutine intersection(p1,p2,p3,p4,q)
    implicit none
    
    real(8), intent(in) :: p1(2), p2(2), p3(2), p4(2)
    real(8), intent(out) :: q
    real(8) :: Det, cross(2)

    print *, 'setting q to -1'
    q = -1._8

    print *, 'calculating D'
    Det = (p1(1)-p2(1))*(p3(2)-p4(2)) - (p1(2)-p2(2))*(p3(1)-p4(1))

    if ( Det == 0 ) then
        print *, 'D is zero'
        return
    end if

    cross(1) = ((p3(1)-p4(1))*(p1(1)*p2(2)-p1(2)*p2(1))-(p1(1)-p2(1))*(p3(1)*p4(2)-p3(2)*p4(1)))/Det
    cross(2) = ((p3(2)-p4(2))*(p1(1)*p2(2)-p1(2)*p2(1))-(p1(2)-p2(2))*(p3(1)*p4(2)-p3(2)*p4(1)))/Det

    !-- check if cross(2) lies inside of the four points pi(2)
    if ( cross(1) >= min(p1(1),p2(1)) .and. cross(1) >= min(p3(1),p4(1)) &
        .and. cross(1) <= max(p1(1),p2(1)) .and. cross(1) <= max(p3(1),p4(1)) &
        .and. cross(2) >= min(p1(2),p2(2)) .and. cross(2) >= min(p3(2),p4(2)) &
        .and. cross(2) <= max(p1(2),p2(2)) .and. cross(2) <= max(p3(2),p4(2))) then
        
        print *, 'intersection is inside'
        !-- calculate q
        q = SQRT((p1(1)-cross(1))**2+(p1(2)-cross(2))**2)/SQRT((p1(1)-p2(1))**2+(p1(2)-p2(2))**2)
    end if
    
end subroutine intersection