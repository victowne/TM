Subroutine right(xo,xi)               ! for rk4, right hand side of mhd equations
        
    USE omp_lib
    use variable
    implicit none
    integer i,j,k,jp1,jm1,ip1,im1
    real*8, dimension(ni,nj,7) :: xi,xo ! 1 2 3 4 5 6 7 --> rho p ux uy uz psi bz
    real*8, dimension(ni,nj) :: rrho                                      ! 1/rho
    real*8 s
    real*8 num4(2)
    integer x0,y0
    real*8 eccd
    va=dsqrt(b0*b0/rho0)
    p0=0.5*beta*(b0**2)
    t0=0.5*beta*va*va
    
    do i=1,ni
        s=(i-(ni+1)/2)*dx/bl
        psie(i)=-bl*log(cosh(s))
        bye(i)=tanh(s)
        pe(i)=-tanh(s*bl/pl)+1
        bze(i)=sqrt(2*(1+.5*(bze0**2-bye(i)**2)-pe(i)))
        rhoe(i)=1.
        ! vys(i)=vy0*(-1./(cosh(s*bl/lv))+1)
        vys(i)=vy0
    enddo

    ! calculate Bx, By
    call calcbxy(xi)

    CALL omp_set_num_threads(numthreads)
    !$omp parallel do private(i,j) 
    do i=1,ni
        do j=1,nj
            rrho(i,j)=1.0/xi(i,j,1)                                 ! rrho=1/rho
            pt(i,j)=xi(i,j,2)+0.5*(bx(i,j)**2+by(i,j)**2+xi(i,j,7)**2)           ! pt=p+b^2/2
        enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(i,j) 
    do j=1,nj
        do i=1,ni
            psi(i,j)=xi(i,j,6)
        end do 
    end do 
    !$omp end parallel do 
    num4=maxloc(psi)
    x0=num4(1)
    y0=num4(2)


    !$omp parallel private(i,j)
    !$omp do
    do j=1,nj
        jp1=j+1
        jm1=j-1
        if(j==nj) jp1=1
        if(j==1) jm1=nj
        do i=2,ni-1
            ip1=i+1
            im1=i-1

            xo(i,j,1)=-xi(i,j,3)*dx2*(xi(ip1,j,1)-xi(im1,j,1))&
                    -xi(i,j,4)*dy2*(xi(i,jp1,1)-xi(i,jm1,1))&
                    -xi(i,j,1)*(dx2*(xi(ip1,j,3)-xi(im1,j,3))&
                    +dy2*(xi(i,jp1,4)-xi(i,jm1,4)))&
                    +drho*(ddx*(xi(ip1,j,1)+xi(im1,j,1)-2.0*xi(i,j,1))&
                    +ddy*(xi(i,jp1,1)+xi(i,jm1,1)-2.0*xi(i,j,1)))   

            xo(i,j,2)=-xi(i,j,3)*dx2*(xi(ip1,j,2)-xi(im1,j,2))&
                    -xi(i,j,4)*dy2*(xi(i,jp1,2)-xi(i,jm1,2))&
                    -gamma*xi(i,j,2)*(dx2*(xi(ip1,j,3)-xi(im1,j,3))&
                    +dy2*(xi(i,jp1,4)-xi(i,jm1,4)))&
                    +dp*(ddx*(xi(ip1,j,2)+xi(im1,j,2)-2.0*xi(i,j,2))&
                    +ddy*(xi(i,jp1,2)+xi(i,jm1,2)-2.0*xi(i,j,2)) &
                    -ddx*(pe(ip1)+pe(im1)-2.*pe(i)))                

            xo(i,j,3)=-xi(i,j,3)*dx2*(xi(ip1,j,3)-xi(im1,j,3))&
                    -xi(i,j,4)*dy2*(xi(i,jp1,3)-xi(i,jm1,3))&
                    +rrho(i,j)*((bx(i,j)*dx2*(bx(ip1,j)-bx(im1,j))&
                    +by(i,j)*dy2*(bx(i,jp1)-bx(i,jm1))&
                    -dx2*(pt(ip1,j)-pt(im1,j)))&
                    +vis*(ddx*(xi(ip1,j,3)+xi(im1,j,3)-2.0*xi(i,j,3))&
                    +ddy*(xi(i,jp1,3)+xi(i,jm1,3)-2.0*xi(i,j,3))))  

            xo(i,j,4)=-xi(i,j,3)*dx2*(xi(ip1,j,4)-xi(im1,j,4))&
                    -xi(i,j,4)*dy2*(xi(i,jp1,4)-xi(i,jm1,4))&
                    +rrho(i,j)*((bx(i,j)*dx2*(by(ip1,j)-by(im1,j))&
                    +by(i,j)*dy2*(by(i,jp1)-by(i,jm1))&
                    -dy2*(pt(i,jp1)-pt(i,jm1)))&
                    +vis*(ddx*(xi(ip1,j,4)+xi(im1,j,4)-2.0*xi(i,j,4))&
                    +ddy*(xi(i,jp1,4)+xi(i,jm1,4)-2.0*xi(i,j,4))))

            xo(i,j,5)=-xi(i,j,3)*dx2*(xi(ip1,j,5)-xi(im1,j,5))&
                    -xi(i,j,4)*dy2*(xi(i,jp1,5)-xi(i,jm1,5))&
                    +rrho(i,j)*((bx(i,j)*dx2*(xi(ip1,j,7)-xi(im1,j,7))&
                    +by(i,j)*dy2*(xi(i,jp1,7)-xi(i,jm1,7)))&                        
                    +vis*(ddx*(xi(ip1,j,5)+xi(im1,j,5)-2.0*xi(i,j,5))&
                    +ddy*(xi(i,jp1,5)+xi(i,jm1,5)-2.0*xi(i,j,5))))

            if(tt.gt.eccdon)then
                    eccd=ejz*exp((-((i-x0)*dx)**2-((j-y0)*dy)**2)/delta2)
            else
                    eccd=0.0
            endif

            xo(i,j,6)=-xi(i,j,3)*dx2*(xi(ip1,j,6)-xi(im1,j,6))&
                    -xi(i,j,4)*dy2*(xi(i,jp1,6)-xi(i,jm1,6))&
                    +eta*(ddx*(xi(ip1,j,6)+xi(im1,j,6)-2.0*xi(i,j,6))&
                    +ddy*(xi(i,jp1,6)+xi(i,jm1,6)-2.0*xi(i,j,6))&
                    +jbs*dx2*(xi(ip1,j,2)-xi(im1,j,2))&
                    -ddx*(psie(ip1)+psie(im1)-2.*psie(i))-eccd)

            xo(i,j,7)=-xi(i,j,3)*dx2*(xi(ip1,j,7)-xi(im1,j,7))&
                    -xi(i,j,4)*dy2*(xi(i,jp1,7)-xi(i,jm1,7))&
                    +bx(i,j)*dx2*(xi(ip1,j,5)-xi(im1,j,5))&
                    +by(i,j)*dy2*(xi(i,jp1,5)-xi(i,jm1,5))&
                    -xi(i,j,7)*(dx2*(xi(ip1,j,3)-xi(im1,j,3))&
                    +dy2*(xi(i,jp1,4)-xi(i,jm1,4)))&                                
                    +eta*(ddx*(xi(ip1,j,7)+xi(im1,j,7)-2.0*xi(i,j,7))&
                    +ddy*(xi(i,jp1,7)+xi(i,jm1,7)-2.0*xi(i,j,7))&
                    -ddx*(bze(ip1)+bze(im1)-2.*bze(i)))             
        enddo
    enddo
    !$omp end do  
    !$omp end parallel 

end Subroutine right
!********************************************************************************
