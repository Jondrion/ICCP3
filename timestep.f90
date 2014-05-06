subroutine timestep(dataarray, x, y, pressure_grad, relaxtime, totaldensity, velocities)
    
    use dispmodule
    integer, intent(in) :: x,y
    real(8), intent(in) :: pressure_grad, relaxtime
    real(8), intent(inout) :: dataarray(y,x,7)
    real(8), intent(out) :: velocities(y,x,2)
    real(8) :: equildensity(y,x,7), totaldensity(y,x)
    integer :: mask(y,x)

    !-- temporary code
    !-- make mask for boundaries: array of all nodes which is 0 for internal and 3 for external points
    mask=0
    mask(1,:)=3
!     mask(:,1)=3
    mask(y,:)=3
!     mask(:,x)=3
    !-- cube in centre
    mask(13:17,10:14)=3

    !-- make density on wall points zero
!     do j = 1, x
!         do i = 1, y
!             if ( mask(i,j) == 3 ) dataarray(i,j,:)=0 
!         end do
!     end do

    call movedensity(dataarray, mask, x, y)

    call calculate_vel(velocities, dataarray, x, y)

    call add_pressure(dataarray,velocities,x,y,pressure_grad)

    call calculate_equildensity(equildensity,totaldensity,velocities, x, y)

!     print *,'equildensity: '
!     call disp(sum(equildensity,3))
    call relax_density(dataarray,equildensity,mask,x,y,relaxtime)


contains

    subroutine movedensity(dataarray,mask,x,y)
        integer, intent(in) :: x,y, mask(y,x)
        real(8), intent(inout) :: dataarray(y,x,7)
        real(8) :: temparray(y,x,7)
        integer :: i,j,k, e_ik(2,7), e_jk(2,7), inew, jnew, knew

        temparray=0

        !-- directions:
        e_ik(1,:)=[0,0,-1,-1,0,1,1]
        e_jk(1,:)=[0,1,1,0,-1,0,1]
        e_ik(2,:)=[0,0,-1,-1,0,1,1]
        e_jk(2,:)=[0,1,0,-1,-1,-1,0]

        do j=1, x
            do i=1,y
                temparray(i,j,1)=dataarray(i,j,1)
                do k=2,7
                    inew=i+e_ik(1+modulo(i,2),k)
                    !-- periodic bc in x-direction
                    jnew=modulo((j+e_jk(1+modulo(i,2),k)-1),x)+1
                    !jnew=j+e_jk(1+modulo(i,2),k)
                    !-- reverse direction if at boundary point
                    knew=modulo((k-2+mask(inew,jnew)),6)+2

                    !-- do bounce-back in one time step and make wall points have zero density
                    if ( mask(inew,jnew) == 3 ) then
                        inew=i
                        jnew=j
                    end if

                    !-- only move densities in direction of domain
                     if ( mask(i,j) == 0 ) then
                        temparray(inew,jnew,knew)=dataarray(i,j,k)
                     end if
                end do
            end do
        end do
        
        dataarray=temparray

    end subroutine movedensity

    subroutine calculate_vel(velocities,dataarray,x,y)
        integer, intent(in) :: x, y
        real(8), intent(in) :: dataarray(y,x,7)
        real(8), intent(out) :: velocities(y,x,2)
        real(8), parameter :: Pi=4*atan(1d0)
        
        integer :: i,j

        do i=1,x
            do j=1,y
                if ( sum(dataarray(j,i,:)) > 0 ) then
                    velocities(j,i,1)=(dataarray(j,i,2)+dataarray(j,i,3)*COS(Pi/3)+dataarray(j,i,4)*COS(2*Pi/3)-dataarray(j,i,5) &
                        -dataarray(j,i,6)*COS(Pi/3)-dataarray(j,i,7)*COS(2*Pi/3))/sum(dataarray(j,i,:))
                    velocities(j,i,2)=(-dataarray(j,i,3)*SIN(Pi/3)-dataarray(j,i,4)*SIN(2*Pi/3) &
                        -dataarray(j,i,6)*SIN(4*Pi/3)-dataarray(j,i,7)*SIN(5*Pi/3))/sum(dataarray(j,i,:))
                end if     
            end do
        end do

    end subroutine calculate_vel

    subroutine add_pressure(dataarray,velocities,x,y,pressure_grad)
        integer, intent(in) :: x,y
        real(8), intent(in) :: pressure_grad, dataarray(y,x,7)
        real(8), intent(inout) :: velocities(y,x,2)
        integer :: i,j

        do i=1,x
            do j=1,y
                if ( sum(dataarray(j,i,:)) > 0 ) velocities(j,i,1)=velocities(j,i,1) + pressure_grad/sum(dataarray(j,i,:))
            end do
        end do

    end subroutine add_pressure

    subroutine relax_density(dataarray,equildensity,mask,x,y,relaxtime)
        integer, intent(in) :: x,y, mask(y,x)
        real(8), intent(in) :: relaxtime, equildensity(y,x,7)
        real(8), intent(inout) :: dataarray(y,x,7)

        !-- relax densities only on internal points (where mask is equal to zero)
        do i = 1, x
            do j = 1, y
                if ( mask(j,i) == 0 ) dataarray(j,i,:) = (1._8-1._8/relaxtime)*dataarray(j,i,:)+(1._8/relaxtime)*equildensity(j,i,:) 
            end do
        end do
!         dataarray(2:y-1,2:x-1,:) = (1._8-1._8/relaxtime)*dataarray(2:y-1,2:x-1,:)+(1._8/relaxtime)*equildensity(2:y-1,2:x-1,:)
        
    end subroutine relax_density


end subroutine timestep

