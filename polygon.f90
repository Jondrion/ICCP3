

subroutine polygon(X_object,np,Object,x,y)
    use dispmodule
    implicit none
    integer, intent(in) :: np,x,y
    real(8), intent(in) :: X_object(np,2)
    real(8) :: X_nodes(y,x,2), Crosspoint(2)
    integer, intent(out) :: Object(y,x)
    integer :: boolean, i, j, k, l, xx, yy, e_ik(2,7), e_jk(2,7), inew, jnew
    real(8) :: q(y,x,6)


    do i=1,x
        xx = i
        do j=1,y
          yy = j
          X_nodes(j,i,1) = xx+modulo(j+1,2)*0.5
          X_nodes(j,i,2) = yy
        end do
    end do
!     print *, "X_object"
!     call disp(X_object)
!     print *, "nodes"
!     call disp(X_nodes(:,:,1))


    Object=0
    q=0
    

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

    Object=(-3)*(Object-1)

!     print *, "Object"
!     call disp(Object)

        e_ik(1,:)=[0,0,-1,-1,0,1,1]
        e_jk(1,:)=[0,1,1,0,-1,0,1]
        e_ik(2,:)=[0,0,-1,-1,0,1,1]
        e_jk(2,:)=[0,1,0,-1,-1,-1,0]

    do j=1, x
        do i=1,y
            if (Object(i,j)==0) then
                do k=2,7
                    inew=i+e_ik(1+modulo(i,2),k)
                    !-- periodic bc in x-direction
                    jnew=modulo((j+e_jk(1+modulo(i,2),k)-1),x)+1
                    if (Object(inew,jnew)==3) then
                        print *, "Check", X_nodes(i,j,:)
                        Crosspoint=0
                        do l=1,np
                            call cross(X_nodes(i,j,:),X_nodes(inew,jnew,:),X_object(l,:),X_object(modulo(l,np)+1,:),Crosspoint)
                            if (e_ik(1+modulo(i,2),k)/=0 .and. Crosspoint(1)/=0) then
                                print *, "eji", e_ik(1+modulo(i,2),k)
                                q(i,j,k-1)=(Crosspoint(1)-X_nodes(i,j,1))/e_ik(1+modulo(i,2),k)
                                exit
                            elseif(Crosspoint(1)/=0) then 
                                print *, "eji2", e_jk(1+modulo(i,2),k)
                                q(i,j,k-1)=(Crosspoint(2)-X_nodes(i,j,2))/e_jk(1+modulo(i,2),k)
                                exit
                            end if
                        end do
                    end if
                end do
            end if
        end do
    end do

    print *, "qvalues"
    call disp( sum(q,3))

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

    subroutine cross(X_node1,X_node2,x1,x2,Crosspoint)
        real(8), intent(in) :: X_node1(2),X_node2(2), x1(2), x2(2)
        real(8), intent(inout) :: Crosspoint(2)
        real(8) :: Diff(2), a1, a2, b1, b2

        Diff=x2-x1
        
        if (Crosspoint(1)/=0) then
            print *, "return"
            return
        end if


        
        if (DIFF(2)==0) then

            a1=(x2(2)-x1(2))/(x2(1)-x1(1))
            a2=(X_node2(2)-X_node1(2))/(X_node2(1)-X_node1(1))

            b1=-(x2(2)-x1(2))/(x2(1)-x1(1))*x1(1)+x1(2)
            b2=-(X_node2(2)-X_node1(2))/(X_node2(1)-X_node1(1))*X_node1(1)+X_node1(2)

            Crosspoint(1)=(b1-b2)/(a2-a1)

            
            if(SIGN(Crosspoint(1),Diff(1))>SIGN(x1(1),Diff(1)) .and. SIGN(Crosspoint(1),Diff(1))<SIGN(x2(1),Diff(1))) then
                Crosspoint(2)=a1*Crosspoint(1)+b1
                print *, "hit1", Crosspoint

            else
                Crosspoint=0               
            end if
        else

            

            a1=(x2(1)-x1(1))/(x2(2)-x1(2))
            b1=-(x2(1)-x1(1))/(x2(2)-x1(2))*x1(2)+x1(1)

            if (X_node1(2)-X_node2(2)==0) then
                b2=X_node2(2)
                Crosspoint(2)=b2
            else
                a2=(X_node2(1)-X_node1(1))/(X_node2(2)-X_node1(2))           
                b2=-(X_node2(1)-X_node1(1))/(X_node2(2)-X_node1(2))*X_node1(2)+X_node1(1)
                Crosspoint(2)=(b1-b2)/(a2-a1)
            end if

            
            if(SIGN(Crosspoint(2),Diff(2))>SIGN(x1(2),Diff(2)) .and. SIGN(Crosspoint(2),Diff(2))<SIGN(x2(2),Diff(2))) then
                Crosspoint(1)=a1*Crosspoint(2)+b1 
                print *, "hit2", Crosspoint               
            else
                Crosspoint=0               
            end if
        end if

    end subroutine cross






end subroutine polygon