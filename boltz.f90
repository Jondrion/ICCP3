program Boltz
    
    use dispmodule

    implicit none
    real(8) :: relaxtime, tubewidth, tubelength, gridsize, pressure, pressure_grad
    real(8), allocatable :: gridarray(:,:,:), velocities(:,:,:), rho(:,:)
    integer :: gridtype, n_y, n_x, i

    call parameters(relaxtime, tubewidth, tubelength, gridsize, gridtype, pressure)
    pressure_grad = pressure/tubelength
    print *, "relaxtime=",relaxtime, "tubewidth", tubewidth, "tubelength", tubelength, "pressure", pressure


!-- Create the grid array
    n_y = int(tubewidth/(gridsize*SQRT(3._8)/2._8))+1
    n_x = int(tubelength/gridsize)+1

    allocate(gridarray(n_y,n_x,gridtype))
    print *,'grid dimensions: ', 'x',SIZE(gridarray,2), 'y:',SIZE(gridarray,1), SIZE(gridarray,3)

    gridarray=0

!-- Calculate intial density distribution (given certain mean velocity)
    allocate(velocities(n_y,n_x,2))
    allocate(rho(n_y,n_x))
    velocities=0._8
    velocities(:,:,1)=0.1_8
    rho=1._8
    call calculate_equildensity(gridarray,rho,velocities,n_x,n_y)
    gridarray(1,:,:)=0
    gridarray(:,1,:)=0
    gridarray(n_y,:,:)=0
    gridarray(:,n_x,:)=0

    do i = 1, 50
      call timestep(gridarray, n_x, n_y, pressure_grad, relaxtime, rho)
      print *,"after timestep ", i, " total density: ", sum(gridarray)
    end do

end program


