Subroutine boundary(xi)   
    ! y direction, periodic boundary condition
    ! x direction boundary condition
    use variable
    implicit none
    integer j
    real*8, dimension(ni,nj,7) :: xi ! 1 2 3 4 5 6 7 --> rho p ux uy uz psi bz

    do j=1,nj         
        xi(1,j,1)=xi(2,j,1)
        xi(1,j,2)=xi(2,j,2)
        xi(1,j,3)=0.
        xi(1,j,4)=xi(2,j,4)
        xi(1,j,5)=xi(2,j,5)
        xi(1,j,6)=2.*xi(2,j,6)-xi(3,j,6)
        ! xi(1,j,6)=0.
        xi(1,j,7)=xi(2,j,7)

        xi(ni,j,1)=xi(ni-1,j,1)
        xi(ni,j,2)=xi(ni-1,j,2)
        ! xi(ni,j,2)=-pper*(1-exp(-tt/tau))*sin(2.*pi*(j-1)*dy/4.)+pp0
        xi(ni,j,3)=0.
        ! xi(ni,j,3)=f0*(exp(-tt/tau)*sin(2.*pi*(j-1)*dy/4.))/tau
        xi(ni,j,4)=xi(ni-1,j,4)
        xi(ni,j,5)=xi(ni-1,j,5)
        xi(ni,j,6)=2.*xi(ni-1,j,6)-xi(ni-2,j,6)
        ! xi(ni,j,6)=psip*sin(2.*pi*(j-1)*dy/4.)*tanh(tt/tau)
        ! if(width*dx.gt.widthon)then
        ! if(tt.gt.widthon)then
        !       xi(ni,j,6)=psiper*sin(mm*2.*(j-1)*dy*pi/4+phi)*tanh((tt-widthon)/ttrmp)!rmp
        !       ! xi(ni,j,6)=-psip*(1-exp(-tt/tau))*sin(mm*2.*pi*(j-1)*dy/4.)
        ! else 
        !       ! xi(ni,j,6)=-psip*(1-exp(-tt/tau))*sin(mm*2.*pi*(j-1)*dy/4.)!force
        !       xi(ni,j,6)=2.*xi(ni-1,j,6)-xi(ni-2,j,6)
        ! end if  
        ! xi(ni,j,6)=psip*tanh(tt/tau)*sin(mm*2.*pi*(j-1)*dy/4.)*(1-tanh((tt-tt0)/tau))
        ! xi(ni,j,6)=psip*tanh(tt/tau)*sin(mm*2.*pi*(j-1)*dy/4.)
        xi(ni,j,7)=xi(ni-1,j,7)
    end do
    return
end Subroutine boundary

!********************************************************************************
