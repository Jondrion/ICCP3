subroutine parameters(lsize, relaxtime, tubewidth, tubelength)

  real(8), intent(out) :: lsize, relaxtime, tubewidth, tubelength
  open(10,file="parameters.txt")
  read(10,*) lsize
  read(10,*) relaxtime
  read(10,*) tubewidth
  read(10,*) tubelength
  close(10)

end subroutine
