program Boltz
    
    use dispmodule

    implicit none
    real(8) :: relaxtime, tubewidth, tubelength, gridsize, pressure, pressure_grad
    real(8), allocatable :: gridarray(:,:,:)
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

    gridarray(3,2,1:7)=1
    gridarray(2,5,1:7)=1

    print *,"intial distribution"
    call disp(sum(gridarray,3))

    do i = 1, 20
      call timestep(gridarray, n_x, n_y, pressure_grad, relaxtime)
      print *,"after timestep ", i, " total density: ", sum(gridarray)
      call disp(sum(gridarray,3))
    end do

end program


