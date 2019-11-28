Subroutine rk4                                ! 4-th Runge-Kutta time integration
    ! required subroutine right(xo,xi)
    !      xi, variable; xo, right side of mhd equations
    ! input: x -- variable array
    ! output: x -- after one time step
    !         kw -- warning, =0 normal, =1 divergent
    USE omp_lib
    use variable
    implicit none
    integer i,j,k
    real*8 s
    
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
        vxs(i)=0.
        ! vys(i)=vy0*cos(0.25*pi*s*bl)
        ! vys(i)=vy0*tanh(s*bl/lv)
        vys(i)=vy0
        ! vys(i)=vy0*(-1./(cosh(s*bl/lv))+1)
        vzs(i)=0.
    enddo
    
    call right(f1,x)
    CALL omp_set_num_threads(numthreads)
    !$omp parallel do private(i,j) 
    do i=2,ni-1
        do j=1,nj       
            y(i,j,1)=x(i,j,1)+0.5*dt*f1(i,j,1)
            y(i,j,2)=x(i,j,2)+0.5*dt*f1(i,j,2)
            y(i,j,3)=x(i,j,3)+0.5*dt*f1(i,j,3)
            y(i,j,4)=x(i,j,4)+0.5*dt*f1(i,j,4)
            y(i,j,5)=x(i,j,5)+0.5*dt*f1(i,j,5)
            y(i,j,6)=x(i,j,6)+0.5*dt*f1(i,j,6)
            y(i,j,7)=x(i,j,7)+0.5*dt*f1(i,j,7)

            y(i,j,1)=y(i,j,1)-rhoe(i)
            y(i,j,2)=y(i,j,2)-pe(i)
            y(i,j,3)=y(i,j,3)-vxs(i)
            y(i,j,4)=y(i,j,4)-vys(i)
            y(i,j,5)=y(i,j,5)-vzs(i)
            y(i,j,6)=y(i,j,6)-psie(i)
            y(i,j,7)=y(i,j,7)-bze(i)
            ! WRITE(*,*) omp_get_thread_num()
        end do
    end do
    !$omp end parallel do
    call boundary(y)

    !$omp parallel do private(i,j) 
    do i=1,ni
        do j=1,nj       
            y(i,j,1)=y(i,j,1)+rhoe(i)
            y(i,j,2)=y(i,j,2)+pe(i)
            y(i,j,3)=y(i,j,3)+vxs(i)
            y(i,j,4)=y(i,j,4)+vys(i)
            y(i,j,5)=y(i,j,5)+vzs(i)                
            y(i,j,6)=y(i,j,6)+psie(i)
            y(i,j,7)=y(i,j,7)+bze(i)
        end do
    end do
    !$omp end parallel do

    call smooth(y)

    call right(f2,y)

    !$omp parallel do private(i,j) 
    do i=2,ni-1
        do j=1,nj       
            y(i,j,1)=x(i,j,1)+0.5*dt*f2(i,j,1)
            y(i,j,2)=x(i,j,2)+0.5*dt*f2(i,j,2)
            y(i,j,3)=x(i,j,3)+0.5*dt*f2(i,j,3)
            y(i,j,4)=x(i,j,4)+0.5*dt*f2(i,j,4)
            y(i,j,5)=x(i,j,5)+0.5*dt*f2(i,j,5)
            y(i,j,6)=x(i,j,6)+0.5*dt*f2(i,j,6)
            y(i,j,7)=x(i,j,7)+0.5*dt*f2(i,j,7)

            y(i,j,1)=y(i,j,1)-rhoe(i)
            y(i,j,2)=y(i,j,2)-pe(i)
            y(i,j,3)=y(i,j,3)-vxs(i)
            y(i,j,4)=y(i,j,4)-vys(i)
            y(i,j,5)=y(i,j,5)-vzs(i)
            y(i,j,6)=y(i,j,6)-psie(i)
            y(i,j,7)=y(i,j,7)-bze(i)
        end do
    end do
    !$omp end parallel do

    call boundary(y)

    !$omp parallel do private(i,j) 
    do i=1,ni
        do j=1,nj       
            y(i,j,1)=y(i,j,1)+rhoe(i)
            y(i,j,2)=y(i,j,2)+pe(i)
            y(i,j,3)=y(i,j,3)+vxs(i)
            y(i,j,4)=y(i,j,4)+vys(i)
            y(i,j,5)=y(i,j,5)+vzs(i)
            y(i,j,6)=y(i,j,6)+psie(i)
            y(i,j,7)=y(i,j,7)+bze(i)
        end do
    end do
    !$omp end parallel do
    call smooth(y)

    call right(f3,y)

    !$omp parallel do private(i,j) 
    do i=2,ni-1
        do j=1,nj       
            y(i,j,1)=x(i,j,1)+dt*f3(i,j,1)
            y(i,j,2)=x(i,j,2)+dt*f3(i,j,2)
            y(i,j,3)=x(i,j,3)+dt*f3(i,j,3)
            y(i,j,4)=x(i,j,4)+dt*f3(i,j,4)
            y(i,j,5)=x(i,j,5)+dt*f3(i,j,5)
            y(i,j,6)=x(i,j,6)+dt*f3(i,j,6)
            y(i,j,7)=x(i,j,7)+dt*f3(i,j,7)

            y(i,j,1)=y(i,j,1)-rhoe(i)
            y(i,j,2)=y(i,j,2)-pe(i)
            y(i,j,3)=y(i,j,3)-vxs(i)
            y(i,j,4)=y(i,j,4)-vys(i)
            y(i,j,5)=y(i,j,5)-vzs(i)
            y(i,j,6)=y(i,j,6)-psie(i)
            y(i,j,7)=y(i,j,7)-bze(i)
        end do
    end do
    !$omp end parallel do

    call boundary(y)

    !$omp parallel do private(i,j) 
    do i=1,ni
        do j=1,nj       
            y(i,j,1)=y(i,j,1)+rhoe(i)
            y(i,j,2)=y(i,j,2)+pe(i)
            y(i,j,3)=y(i,j,3)+vxs(i)
            y(i,j,4)=y(i,j,4)+vys(i)
            y(i,j,5)=y(i,j,5)+vzs(i)
            y(i,j,6)=y(i,j,6)+psie(i)
            y(i,j,7)=y(i,j,7)+bze(i)
        end do
    end do
    !$omp end parallel do

    call smooth(y)

    call right(f4,y)

    !$omp parallel do private(i,j)
    do i=2,ni-1
        do j=1,nj       
            x(i,j,1)=x(i,j,1)+dt*(f1(i,j,1)+2.*f2(i,j,1)+2.*f3(i,j,1)+f4(i,j,1))/6.
            x(i,j,2)=x(i,j,2)+dt*(f1(i,j,2)+2.*f2(i,j,2)+2.*f3(i,j,2)+f4(i,j,2))/6.
            x(i,j,3)=x(i,j,3)+dt*(f1(i,j,3)+2.*f2(i,j,3)+2.*f3(i,j,3)+f4(i,j,3))/6.
            x(i,j,4)=x(i,j,4)+dt*(f1(i,j,4)+2.*f2(i,j,4)+2.*f3(i,j,4)+f4(i,j,4))/6.
            x(i,j,5)=x(i,j,5)+dt*(f1(i,j,5)+2.*f2(i,j,5)+2.*f3(i,j,5)+f4(i,j,5))/6.
            x(i,j,6)=x(i,j,6)+dt*(f1(i,j,6)+2.*f2(i,j,6)+2.*f3(i,j,6)+f4(i,j,6))/6.
            x(i,j,7)=x(i,j,7)+dt*(f1(i,j,7)+2.*f2(i,j,7)+2.*f3(i,j,7)+f4(i,j,7))/6.

            x(i,j,1)=x(i,j,1)-rhoe(i)
            x(i,j,2)=x(i,j,2)-pe(i)
            x(i,j,3)=x(i,j,3)-vxs(i)
            x(i,j,4)=x(i,j,4)-vys(i)
            x(i,j,5)=x(i,j,5)-vzs(i)
            x(i,j,6)=x(i,j,6)-psie(i)
            x(i,j,7)=x(i,j,7)-bze(i)
        end do
    end do
    !$omp end parallel do

    call boundary(x)

    !$omp parallel do private(i,j) 
    do i=1,ni
        do j=1,nj       
            x(i,j,1)=x(i,j,1)+rhoe(i)
            x(i,j,2)=x(i,j,2)+pe(i)
            x(i,j,3)=x(i,j,3)+vxs(i)
            x(i,j,4)=x(i,j,4)+vys(i)
            x(i,j,5)=x(i,j,5)+vzs(i)
            x(i,j,6)=x(i,j,6)+psie(i)
            x(i,j,7)=x(i,j,7)+bze(i)
        end do
    end do
    !$omp end parallel do 

    call smooth(x)
    
    do k=1,7                                             ! judge divergent or not
        do j=1,nj
            do i=1,ni
                if(abs(x(i,j,k))>large) then
                    kw=1
                    return
                endif
            enddo
        enddo
    enddo
    
    return
end Subroutine rk4

!********************************************************************************
