subroutine parameters(relaxtime, tubewidth, tubelength, gridsize, gridtype, pressure)

  real(8), intent(out) :: relaxtime, tubewidth, tubelength, gridsize, pressure
  integer, intent(out) :: gridtype
  open(10,file="parameters.txt")
  read(10,*) relaxtime
  read(10,*) tubewidth
  read(10,*) tubelength
  read(10,*) gridsize
  read(10,*) gridtype
  read(10,*) pressure
  close(10)

end subroutine
