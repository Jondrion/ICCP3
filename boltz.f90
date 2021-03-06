program Boltz
    
    use dispmodule
    use PLplot3D

    implicit none
    real(8) :: relaxtime, tubewidth, tubelength, gridsize, pressure, pressure_grad, V_object(2), M_object, &
        CoM(2), I_object, alpha_object, object_size, offset(2)
    real(8), allocatable :: gridarray(:,:,:), velocities(:,:,:), rho(:,:), X_object(:,:)
    integer :: gridtype, n_y, n_x, i, n_vertices
    
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

!-- Initialise X_object
    n_vertices=4
    allocate(X_object(n_vertices,2))

!-- Calculate intial density distribution (given certain mean velocity)
    allocate(velocities(n_y,n_x,2))
    allocate(rho(n_y,n_x))

    velocities=0._8
    velocities(:,:,1)=0.1_8
    rho=1._8
    call calculate_equildensity(gridarray,rho,velocities,n_x,n_y)
!     print *,'initial distribution: '
!     call disp(sum(gridarray,3))


    !X_object(1,:)=[1.5_8,3._8/2._8*sqrt(3._8)+5.5_8]
    !X_object(2,:)=[3.5_8,3._8/2._8*sqrt(3._8)+5.5_8]
    !X_object(3,:)=[3.5_8,3._8/2._8*sqrt(3._8)+10.5_8]
    !X_object(4,:)=[1.5_8,3._8/2._8*sqrt(3._8)+10.5_8]
!     X_object(1,:)=[0._8,0._8]
!     X_object(2,:)=[0._8,0._8]
!     X_object(3,:)=[0._8,0._8]
!     X_object(4,:)=[0._8,0._8]
    V_object(1)=0
    V_object(2)=0
    M_object=200
    I_object=3000

    CoM(1)=10.5_8
    CoM(2)=3._8/2._8*sqrt(3._8)+13.0_8

    alpha_object=0
    object_size=3._8
    offset=[0._8,0._8]

    call make_object(n_vertices, CoM, object_size, offset, X_object)


    do i = 1, 400
      call timestep(gridarray, n_x, n_y, pressure_grad, relaxtime, rho, &
        X_object, n_vertices, V_object, M_object, I_object, CoM,alpha_object, velocities)
      print *,"after timestep ", i, " total density: ", sum(gridarray), &
        "velox", minval(velocities(:,:,1)), "veloy", minval(velocities(:,:,2))
      call plot_points(velocities, n_x, n_y, X_object, n_vertices, CoM)
    end do
    
     
end program


