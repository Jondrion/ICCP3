subroutine timestep(dataarray, x, y, pressure_grad, relaxtime)
    
    use dispmodule
    integer, intent(in) :: x,y
    real(8), intent(in) :: pressure_grad, relaxtime
    real(8), intent(inout) :: dataarray(y,x,7)
    real(8) :: velocities(y,x,2)
    real(8) :: equildensity(y,x,7)

    

    call movedensity(dataarray,x,y)
    
    call reverse_bnd_vel(dataarray,x,y)

    call calculate_vel(velocities, dataarray, x, y)

    call add_pressure(dataarray,velocities,x,y,pressure_grad)

    call calculate_equildensity(equildensity,dataarray,velocities, x, y)

    call relax_density(dataarray,equildensity,x,y,relaxtime)


contains

    subroutine movedensity(dataarray, x, y)
        integer, intent(in) :: x,y
        real(8), intent(inout) :: dataarray(y,x,7)
        real(8) :: temparray(y,x,7)
        integer :: i,j,k

        temparray=0

        do j=1, x
            do i=1,y
                do k=1,7
                    if (dataarray(i,j,k)==0) then
                        !- do nothing, you are probably at a boundary
                    else if (mod(i,2)==0) then
                        if (k==1) then
                            temparray(i,j,k)=dataarray(i,j,k)
                        else if (k==2) then
                            temparray(i,j+1,k)=dataarray(i,j,k)
                        else if (k==3) then
                            temparray(i-1,j+1,k)=dataarray(i,j,k)
                        else if (k==4) then
                            temparray(i-1,j,k)=dataarray(i,j,k)
                        else if (k==5) then
                            temparray(i,j-1,k)=dataarray(i,j,k)
                        else if (k==6) then
                            temparray(i+1,j,k)=dataarray(i,j,k)
                        else if (k==7) then
                            temparray(i+1,j+1,k)=dataarray(i,j,k)
                        end if
                    else
                        if (k==1) then
                            temparray(i,j,k)=dataarray(i,j,k)
                        else if (k==2) then
                            temparray(i,j+1,k)=dataarray(i,j,k)
                        else if (k==3) then
                            temparray(i-1,j,k)=dataarray(i,j,k)
                        else if (k==4) then
                            temparray(i-1,j-1,k)=dataarray(i,j,k)
                        else if (k==5) then
                            temparray(i,j-1,k)=dataarray(i,j,k)
                        else if (k==6) then
                            temparray(i+1,j-1,k)=dataarray(i,j,k)
                        else if (k==7) then
                            temparray(i+1,j,k)=dataarray(i,j,k)
                        end if
                    end if
                end do
            end do
        end do
        
        dataarray=temparray

    end subroutine movedensity

    subroutine reverse_bnd_vel(dataarray, x, y)
        integer, intent(in) :: x, y
        real(8), intent(inout) :: dataarray(y,x,7)
        real(8) :: tempval
        
        integer :: i,k

        do i=1,x
            do k=2,4
                !- top row beyond system boundaries
                tempval=dataarray(1,i,k)
                dataarray(1,i,k)=dataarray(1,i,k+3)
                dataarray(1,i,k+3)=tempval

                !- bottom row beyond system boundaries
                tempval=dataarray(y,i,k)
                dataarray(y,i,k)=dataarray(y,i,k+3)
                dataarray(y,i,k+3)=tempval
            end do
        end do

        do i=2,y-1
            do k=2,4
                !- left row beyond system boundaries
                tempval=dataarray(i,1,k)
                dataarray(i,1,k)=dataarray(i,1,k+3)
                dataarray(i,1,k+3)=tempval

                !- right row beyond system boundaries
                tempval=dataarray(i,x,k)
                dataarray(i,x,k)=dataarray(i,x,k+3)
                dataarray(i,x,k+3)=tempval
            end do
        end do

    end subroutine reverse_bnd_vel

    subroutine calculate_vel(velocities,dataarray,x,y)
        integer, intent(in) :: x, y
        real(8), intent(in) :: dataarray(y,x,7)
        real(8), intent(out) :: velocities(y,x,2)
        real(8), parameter :: Pi=4*atan(1d0)
        
        integer :: i,j

        do i=1,x
            do j=1,y
                velocities(j,i,1)=(dataarray(j,i,2)+dataarray(j,i,3)*COS(Pi/3)+dataarray(j,i,4)*COS(2*Pi/3)-dataarray(j,i,5) &
                    -dataarray(j,i,6)*COS(Pi/3)-dataarray(j,i,7)*COS(2*Pi/3))
                velocities(j,i,2)=(-dataarray(j,i,3)*SIN(Pi/3)-dataarray(j,i,4)*SIN(2*Pi/3) &
                    -dataarray(j,i,6)*SIN(4*Pi/3)-dataarray(j,i,7)*SIN(5*Pi/3))              
            end do
        end do

    end subroutine calculate_vel

    subroutine add_pressure(dataarray,velocities,x,y,pressure_grad)
        integer, intent(in) :: x,y
        real(8), intent(in) :: pressure_grad, dataarray(y,x,7)
        real(8), intent(inout) :: velocities(y,x,2)
        integer :: i,j

        do i=2,x-1
            do j=2,y-1
                if ( sum(dataarray(j,i,:)) > 0 ) velocities(j,i,1)=velocities(j,i,1) + pressure_grad/sum(dataarray(j,i,:))
            end do
        end do

    end subroutine add_pressure

    subroutine calculate_equildensity(equildensity, dataarray, velocities,x,y)

        use dispmodule
        integer, intent(in) :: x,y
        real(8),  intent(in) :: velocities(y,x,2), dataarray(y,x, 7)
        real(8),  intent(out) :: equildensity(y,x,7)
        real(8) :: direction1D(12)
        real(8) :: direction(6,2)
        real(8), parameter :: Pi=4*atan(1d0)
        real(8) :: totaldensity(y,x)
        integer :: i,j,k

        direction1D=[1._8,COS(Pi/3),COS(2*Pi/3),-1._8,COS(4*Pi/3),COS(5*Pi/3), &
             0._8,SIN(Pi/3),SIN(2*Pi/3),0._8,SIN(4*Pi/3),SIN(5*Pi/3)]
        direction=reshape(direction1D,[6,2])

        totaldensity=sum(dataarray,3)
        equildensity=0

        print *,"totaldensity: "
        call disp(totaldensity)

         do i=1,x
            do j=1,y
                do k=1,7
                    if (k==1) then
                        equildensity(j,i,k)=1._8/2*totaldensity(j,i)*(1-2*dot_product(velocities(j,i,1:2),velocities(j,i,1:2)))                        
                    else
                        equildensity(j,i,k)=1._8/12*totaldensity(j,i)*(1+4*dot_product(velocities(j,i,1:2),direction(k-1,1:2)) &
                            -2*dot_product(velocities(j,i,1:2),velocities(j,i,1:2)) &
                                +8*dot_product(velocities(j,i,1:2),direction(k-1,1:2))**2)
                    end if
                end do
            end do
        end do

        print *,"equildensity: "
        call disp(equildensity(:,:,2))
    
    end subroutine calculate_equildensity

    subroutine relax_density(dataarray,equildensity,x,y,relaxtime)
        integer, intent(in) :: x,y
        real(8), intent(in) :: relaxtime, equildensity(y,x,7)
        real(8), intent(inout) :: dataarray(y,x,7)

        dataarray(2:y-1,2:x-1,:) = (1._8-1._8/relaxtime)*dataarray(2:y-1,2:x-1,:)+(1._8/relaxtime)*equildensity(2:y-1,2:x-1,:)
        
    end subroutine relax_density


end subroutine timestep

