subroutine timestep(dataarray, x, y)
	integer, intent(in) :: x,y
	real(8), intent(inout) :: dataarray(x,y,7)

end subroutine

contains

subroutine movedensity(dataarray, x, y)
	real(8), intent(inout) :: dataarray(x,y,7)
	real(8) :: temparray(x,y,7)
	integer :: i,j,k

	do i=1, x
		do j=1,y
			do k=1,7
				if (k==1) then
					temparray(i,j,k)=temparray(i,j,k)+dataarray(i,j,k)
				else if (k==2) then
					temparray(i,j,k)=temparray(i,j,k)+dataarray(i,j,k)
				else if (k==3) then
					temparray(i,j,k)=temparray(i,j,k)+dataarray(i,j,k)
				else if (k==4) then
					temparray(i,j,k)=temparray(i,j,k)+dataarray(i,j,k)
				else if (k==5) then
					temparray(i,j,k)=temparray(i,j,k)+dataarray(i,j,k)
				else if (k==6) then
					temparray(i,j,k)=temparray(i,j,k)+dataarray(i,j,k)
				end if
			end do
		end do
	end do

end subroutine
