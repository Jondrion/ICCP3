program Boltz

  implicit none
  real(8) :: lsize, relaxtime, tubewidth, tubelength, gridsize
  real(8), allocatable :: gridarray(:,:,:)
  integer :: gridtype, n_y, n_x

  call parameters(lsize, relaxtime, tubewidth, tubelength, gridsize, gridtype)
  print *, "lsize=",lsize, "relaxtime=",relaxtime, "tubewidth", tubewidth, "tubelength", tubelength


!-- Create the grid array
  n_y = int(tubewidth/(gridsize*SQRT(3._8)/2._8))+1
  n_x = int(tubelength/gridsize)+1

  allocate(gridarray(n_y,n_x,gridtype))
  print *,'grid dimensions: ', SIZE(gridarray,1), SIZE(gridarray,2), SIZE(gridarray,3)
  
 ! call timestep(gridarray, n_x, n_y)

end program


