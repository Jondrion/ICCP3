subroutine timestep(dataarray, x, y, pressure_grad, relaxtime, totaldensity, X_object, n_vertices, V_object, &
    M_object, CoM, alpha_object, velocities)
    
    use dispmodule
    implicit none
    integer, intent(in) :: x,y, n_vertices
    real(8), intent(in) :: pressure_grad, relaxtime, M_object
    real(8), intent(inout) :: dataarray(y,x,7), V_object(2), CoM(2)
    real(8), intent(out) :: velocities(y,x,2)
    real(8) :: equildensity(y,x,7), totaldensity(y,x), q(y,x,7)
    real(8), intent(inout) :: X_object(n_vertices,2), alpha_object
    integer :: mask(y,x), Object(y,x)
    real(8) :: DV(y,x,2)
    

    !-- temporary code
    !-- make mask for boundaries: array of all nodes which is 0 for internal and 3 for external points
    mask=0
    mask(1,:)=3
!     mask(:,1)=3
    mask(y,:)=3
!     mask(:,x)=3
    !-- cube in centre
!     mask(14:15,10:12)=3
!     mask(7:27,20:22)=3
 

    call polygon(X_object,4,Object,x,y,q)

    mask(2:y-1,:)=Object(2:y-1,:)

    call movedensity(dataarray, DV, mask, x, y, q)

    call calculate_vel(velocities, dataarray, x, y)

    call add_pressure(dataarray,velocities,x,y,pressure_grad)

    call calculate_equildensity(equildensity,totaldensity,velocities, x, y)

    call relax_density(dataarray,equildensity,mask,x,y,relaxtime)

    
    call moveobject(X_object, n_vertices,V_object,M_object, DV,x,y,CoM, alpha_object)

    




contains

    subroutine movedensity(dataarray, DV, mask,x,y,q)
        integer, intent(in) :: x,y, mask(y,x)
        real(8), intent(inout) :: dataarray(y,x,7)
        real(8), intent(inout) :: q(y,x,7)
        real(8) :: temparray(y,x,7)
        integer :: i,j,k, e_ik(2,7), e_jk(2,7), inew, jnew, knew, i1, i2, j1, j2
        real(8) :: V1(y,x,2), v2(2)
        real(8), intent(out) :: DV(y,x,2)

        temparray=0
        DV=0

        !-- directions:
        e_ik(1,:)=[0,0,-1,-1,0,1,1]
        e_jk(1,:)=[0,1,1,0,-1,0,1]
        e_ik(2,:)=[0,0,-1,-1,0,1,1]
        e_jk(2,:)=[0,1,0,-1,-1,-1,0]

        call calculate_vel(V1, dataarray, x, y)

        do j=1, x
            do i=1,y
                temparray(i,j,1)=dataarray(i,j,1)
                do k=2,7
                    inew=i+e_ik(1+modulo(i,2),k)
                    !-- periodic bc in x-direction
                    jnew=modulo((j+e_jk(1+modulo(i,2),k)-1),x)+1
                    !-- reverse direction if at boundary point
                    knew=modulo((k-2+mask(inew,jnew)),6)+2

                    !-- do bounce-back in one time step and use interpolation
                    if ( mask(inew,jnew) == 3 ) then
                        inew=i
                        jnew=j
                    end if

                    !-- only move densities in direction of domain
                     if ( mask(i,j) == 0 .and. mask(inew,jnew) == 0 ) then
                        temparray(inew,jnew,knew)=dataarray(i,j,k)
                     end if
                end do
            end do
        end do

       
        do j=1,x
            do i=1,y
                do k=2,7
                    inew=i+e_ik(1+modulo(i,2),k)
                    !-- periodic bc in x-direction
                    jnew=modulo((j+e_jk(1+modulo(i,2),k)-1),x)+1
                    !-- reverse direction if at boundary point
                    knew=modulo((k-2+mask(inew,jnew)),6)+2

                    !-- do bounce-back in one time step and use interpolation
                    if ( mask(inew,jnew) == 3 .and. mask(i,j)==0 ) then
                        inew=i
                        jnew=j
                        !-- calculate location of two interpolation neighbours
                        i1=i+e_ik(1+modulo(i,2),knew)
                        j1=modulo((j+e_jk(1+modulo(i,2),knew)-1),x)+1
                        i2=i1+e_ik(1+modulo(i1,2),knew)
                        j2=modulo((j1+e_jk(1+modulo(i1,2),knew)-1),x)+1
                        !-- different interpolation depending on q

                        if ( q(i,j,k) < 0.5_8 ) then

                            temparray(inew,jnew,knew) = ( q(i,j,k)*(1._8+2._8*q(i,j,k))*dataarray(i,j,k) & 
                            + (1._8-4._8*q(i,j,k)**2)*dataarray(i1,j1,k) - q(i,j,k)*(1._8-2._8*q(i,j,k))*dataarray(i2,j2,k) )
                            
                            call calculate_vel(v2,temparray(i,j,:),1,1)
                            DV(i,j,:)=V1(i,j,:)-v2

                        end if

                        if ( q(i,j,k) >= 0.5_8 ) then

                            temparray(inew,jnew,knew) = ( 1._8/(q(i,j,k)*(2*q(i,j,k)+1))*dataarray(i,j,k) &
                            + (2._8*q(i,j,k)-1)/q(i,j,k)*temparray(i1,j1,knew) - (2*q(i,j,k)-1._8)/(2*q(i,j,k)+1) & 
                            *temparray(i2,j2,knew) )
                            
                            call calculate_vel(v2,temparray(i,j,:),1,1)
                            DV(i,j,:)=V1(i,j,:)-v2

                        end if
                    else
                   
                    end if                    
                end do
            end do
        end do
        print *,"Q"
        call disp(sum(q,3))


        !Print *, "DVx"
        !call disp(DV(:,:,1))
        !Print *, "DVy"
        !call disp(DV(:,:,2))
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
        integer :: i,j

        !-- relax densities only on internal points (where mask is equal to zero)
        do i = 1, x
            do j = 1, y
                if ( mask(j,i) == 0 ) dataarray(j,i,:) = (1._8-1._8/relaxtime)*dataarray(j,i,:)+(1._8/relaxtime)*equildensity(j,i,:) 
            end do
        end do
        
    end subroutine relax_density

    subroutine moveobject(X_object, n_vertices, V_object, M_object, DV,x,y,CoM, alpha_object)
        integer, intent(in) :: n_vertices, x, y
        real(8), intent(in) :: DV(y,x,2), M_object
        real(8), intent(inout) :: X_object(n_vertices,2), V_object(2),CoM(2), alpha_object
        real(8) :: X_nodes(y,x,2), Ftotal(2),R(2), I_object, Torq
        integer :: i,j

            Ftotal=sum(sum(DV,1),1)

         

            V_object=V_object+Ftotal/M_object

            X_object(:,1)=X_object(:,1)+V_object(1)
            X_object(:,2)=X_object(:,2)+V_object(2)
            CoM=CoM+V_object

            R=0
            Torq=0
        do i=1,x            
            do j=1,y
                if (DV(j,i,1)/=0 .or. DV(j,i,2)/=0) then
                    X_nodes(j,i,1) = i+modulo(j+1,2)*0.5
                    X_nodes(j,i,2) = j*(sqrt(3._8)/2)

                    R=X_nodes(j,i,:)-CoM
                    Torq=Torq+DV(j,i,2)*R(1)-DV(j,i,1)*R(2)
                end if              
            end do
        end do
        I_object=100
        alpha_object=alpha_object+Torq/I_object

       X_object(:,1)=(X_object(:,1)-CoM(1))*cos(alpha_object)-(X_object(:,2)-CoM(2))*sin(alpha_object)+CoM(1)
       X_object(:,2)=(X_object(:,1)-CoM(1))*sin(alpha_object)+(X_object(:,2)-CoM(2))*cos(alpha_object)+CoM(2)

    end subroutine moveobject

        


end subroutine timestep

