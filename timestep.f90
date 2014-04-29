subroutine timestep(dataarray, x, y)
    integer, intent(in) :: x,y
    real(8), intent(inout) :: dataarray(y,x,7)

    call movedensity(dataarray,x,y)


contains

    subroutine movedensity(dataarray, x, y)
        integer, intent(in) :: x,y
        real(8), intent(inout) :: dataarray(y,x,7)
        real(8) :: temparray(y,x,7)
        integer :: i,j,k

        do i=2, x-1
            do j=2,y-1
                do k=1,7
                    if (mod(i,2)==0) then
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

end subroutine timestep

