subroutine calculate_equildensity(equildensity, totaldensity, velocities,x,y)

    use dispmodule
    integer, intent(in) :: x,y
    real(8),  intent(in) :: velocities(y,x,2)
    real(8),  intent(out) :: equildensity(y,x,7)
    real(8) :: direction1D(12)
    real(8) :: direction(6,2)
    real(8), parameter :: Pi=4*atan(1d0)
    real(8), intent(in) :: totaldensity(y,x)
    integer :: i,j,k

    direction1D=[1._8,COS(Pi/3),COS(2*Pi/3),-1._8,COS(4*Pi/3),COS(5*Pi/3), &
         0._8,SIN(Pi/3),SIN(2*Pi/3),0._8,SIN(4*Pi/3),SIN(5*Pi/3)]
    direction=reshape(direction1D,[6,2])

    equildensity=0

!     print *,"totaldensity: "
!     call disp(totaldensity)

     do i=2,x-1
        do j=2,y-1
            do k=1,7
                if (k==1) then
                    equildensity(j,i,k)=1._8/2*totaldensity(j,i)*(1-2*dot_product(velocities(j,i,1:2),velocities(j,i,1:2)))                        
                else
                    equildensity(j,i,k)=1._8/12*totaldensity(j,i)*(1+4*dot_product(velocities(j,i,1:2),direction(k-1,1:2)) &
                        -2*dot_product(velocities(j,i,1:2),velocities(j,i,1:2)) &
                            +8*dot_product(velocities(j,i,1:2),direction(k-1,1:2))**2)
                end if
            end do
        end do
    end do

!     print *,"equildensity: "
!     call disp(equildensity(:,:,2))

end subroutine calculate_equildensity