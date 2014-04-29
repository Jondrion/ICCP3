subroutine timestep(dataarray, x, y)
	integer, intent(in) :: x,y
	real(8), intent(inout) :: dataarray(x,y,7)

end subroutine

contains

subroutine movedensity(dataarray, x, y)
	real(8), intent(inout) :: dataarray(x,y,7)
	real(8) :: temparray(x,y,7)
	integer :: i,j,k

end subroutine
