program Boltz

  implicit none
  real(8) :: lsize, relaxtime, tubewidth, tubelength

  call parameters(lsize, relaxtime, tubewidth, tubelength)
  print *, "lsize=",lsize, "relaxtime=",relaxtime, "tubewidth", tubewidth, "tubelength", tubelength

  
end program


