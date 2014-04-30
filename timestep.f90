subroutine timestep(dataarray, x, y)
    integer, intent(in) :: x,y
    real(8), intent(inout) :: dataarray(y,x,7)
    real(8) velocities(y,x,2)

    call movedensity(dataarray,x,y)

    call reverse_bnd_vel(dataarray,x,y)


contains

    subroutine movedensity(dataarray, x, y)
        real(8), intent(inout) :: dataarray(y,x,7)
        real(8) :: temparray(y,x,7)
        integer :: i,j,k

        do i=2, x-1
            do j=2,y-1
                do k=1,7
                    if (mod(j,2)==0)
                        if (k==1) then
                            temparray(i,j,k)=dataarray(i,j,k)
                        else if (k==2) then
                            temparray(i+1,j,k)=dataarray(i,j,k)
                        else if (k==3) then
                            temparray(i+1,j-1,k)=dataarray(i,j,k)
                        else if (k==4) then
                            temparray(i,j-1,k)=dataarray(i,j,k)
                        else if (k==5) then
                            temparray(i-1,j,k)=dataarray(i,j,k)
                        else if (k==6) then
                            temparray(i,j+1,k)=dataarray(i,j,k)
                        else if (k==7) then
                            temparray(i+1,j+1,k)=dataarray(i,j,k)
                        end if
                    else
                        if (k==1) then
                            temparray(i,j,k)=dataarray(i,j,k)
                        else if (k==2) then
                            temparray(i+1,j,k)=dataarray(i,j,k)
                        else if (k==3) then
                            temparray(i,j-1,k)=dataarray(i,j,k)
                        else if (k==4) then
                            temparray(i-1,j-1,k)=dataarray(i,j,k)
                        else if (k==5) then
                            temparray(i-1,j,k)=dataarray(i,j,k)
                        else if (k==6) then
                            temparray(i-1,j+1,k)=dataarray(i,j,k)
                        else if (k==7) then
                            temparray(i,j+1,k)=dataarray(i,j,k)
                        end if
                    end if
                end do
            end do
        end do
        
        dataarray=temparray

    end subroutine movedensity

    subroutine reverse_bnd_vel(dataarray, x, y)
        real(8), intent(inout) :: dataarray(y,x,7)
        real(8) :: tempval
        integer, intent(in) :: x, y
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

    end subroutine

    subroutine calculate_vel(velocities,dataarray,x,y)
        real(8), intent(inout) :: dataarray(y,x,7)
        real(8), intent(out) :: velocities(y,x,2)
        integer, intent(in) :: x, y
        integer :: i,j

        do i=1,x
            do j=1,y
                velocities(y,x,1)=(dataarray(y,x,2)+dataarray(y,x,3)*COSD(60)+dataarray(y,x,4)*COSD(120)-dataarray(y,x,5)-dataarray(y,x,6)*COSD(60)-dataarray(y,x,7)*COSD(120))/(sum(dataarray, dim=3))
                velocities(y,x,2)=(dataarray(y,x,2)+dataarray(y,x,3)*SIND(60)+dataarray(y,x,4)*SIND(120)-dataarray(y,x,5)-dataarray(y,x,6)*SIND(60)-dataarray(y,x,7)*SIND(120))/(sum(dataarray, dim=3))
            end do
        end do

end subroutine timestep
