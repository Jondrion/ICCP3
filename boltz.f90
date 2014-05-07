program Boltz
    
    use dispmodule
    use PLplot3D

    implicit none
    real(8) :: relaxtime, tubewidth, tubelength, gridsize, pressure, pressure_grad
    real(8), allocatable :: gridarray(:,:,:), velocities(:,:,:), rho(:,:), Object(:,:)
    integer :: gridtype, n_y, n_x, i
    real(8) :: X_object(4,2)

    call parameters(relaxtime, tubewidth, tubelength, gridsize, gridtype, pressure)
    pressure_grad = pressure/tubelength
    print *, "relaxtime=",relaxtime, "tubewidth", tubewidth, "tubelength", tubelength, "pressure", pressure

    call plot_init()

!-- Create the grid array
    n_y = int(tubewidth/(gridsize*SQRT(3._8)/2._8))+1
    n_x = int(tubelength/gridsize)+1

    allocate(gridarray(n_y,n_x,gridtype))
    print *,'grid dimensions: ', 'x',SIZE(gridarray,2), 'y:',SIZE(gridarray,1), SIZE(gridarray,3)

    gridarray=0

!-- Calculate intial density distribution (given certain mean velocity)
    allocate(velocities(n_y,n_x,2))
    allocate(rho(n_y,n_x))
!     allocate(Object(n_y,n_x))
    velocities=0._8
    velocities(:,:,1)=0.1_8
    rho=1._8
    call calculate_equildensity(gridarray,rho,velocities,n_x,n_y)
!     print *,'initial distribution: '
!     call disp(sum(gridarray,3))

!     X_object(1,:)=[2._8,2._8]
!     X_object(2,:)=[8._8,2._8]
!     X_object(3,:)=[8._8,8._8]
!     X_object(4,:)=[2._8,8._8]
!     call polygon(X_object,4,Object,n_x,n_y)

    do i = 1, 100
      call timestep(gridarray, n_x, n_y, pressure_grad, relaxtime, rho, velocities)
      print *,"after timestep ", i, " total density: ", sum(gridarray)
      call plot_points(velocities, n_x, n_y)
    end do
    
!     print *,'velocity x: '
!     call disp(velocities(:,:,1))
!     print *,'velocity y: '
!     call disp(velocities(:,:,2))

end program


