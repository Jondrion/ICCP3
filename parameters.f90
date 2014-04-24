subroutine parameters(lsize, relaxtime, tubewidth, tubelength, gridsize, gridtype)

  real(8), intent(out) :: lsize, relaxtime, tubewidth, tubelength, gridsize
  integer, intent(out) :: gridtype
  open(10,file="parameters.txt")
  read(10,*) lsize
  read(10,*) relaxtime
  read(10,*) tubewidth
  read(10,*) tubelength
  read(10,*) gridsize
  read(10,*) gridtype
  close(10)

end subroutine
