subroutine polygon(X_object,np,Object,x,y)
    use dispmodule
    integer, intent(in) :: np,x,y
    real(8), intent(in) :: X_object(np,2)
    real(8) :: X_nodes(y,x,2)
    integer, intent(out) :: Object(y,x)
    integer :: boolean, i, j, k, xx, yy


    do i=1,x
        xx = i
        do j=1,y
          yy = j
          X_nodes(j,i,1) = xx+modulo(j+1,2)*0.5
          X_nodes(j,i,2) = yy
        end do
    end do
    print *, "X_object"
    call disp(X_object)
    print *, "nodes"
    call disp(X_nodes(:,:,1))


    Object=0
    

    do i=1,x
        do j=1,y
            boolean=0
            do k=1,np                
                call line(X_object(k,:),X_object(modulo(k,np)+1,:),X_nodes(j,i,:),boolean)
            end do
            !print *,"bool",boolean,"xy",X_nodes(j,i,:)
            Object(j,i)=boolean
        end do
    end do

    print *, "Object"
    call disp(Object)





contains

    subroutine line(x1,x2,X,boolean)
        use dispmodule
        real(8), intent(in) :: x1(2),x2(2),X(2)
        integer, intent(inout) :: boolean
        real(8) :: functionvalue, Diff(2)

        Diff=x2-x1

        if (boolean/=0) then
            return
        end if

        
        
        if (DIFF(2)==0) then
            functionvalue=(x2(2)-x1(2))/(x2(1)-x1(1))*(X(1)-x1(1))+x1(2)
            if(SIGN(X(2),Diff(1))>SIGN(functionvalue,Diff(1))) then
                boolean=0

            else
                boolean=1                
            end if
        else
            functionvalue=(x2(1)-x1(1))/(x2(2)-x1(2))*(X(2)-x1(2))+x1(1)
            if(SIGN(X(1),Diff(2))<SIGN(functionvalue,Diff(2))) then
                boolean=0
            else
                boolean=1              
            end if
        end if

        

    end subroutine line



end subroutine polygon