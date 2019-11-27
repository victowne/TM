Subroutine initial

    use variable
    implicit none
    integer i,j
    real*8 s,sii,sij
    
    va=dsqrt(b0*b0/rho0)
    p0=0.5*beta*(b0*b0)
    t0=0.5*beta*va*va
    ta=bl/va
    tau=numta*ta
    tt0=numta1*ta
    ttrmp=numta2*ta

    x=0.0
    
    do i=1,ni
        s=(i-(ni+1)/2)*dx/bl
        psie(i)=-bl*log(cosh(s))
        bye(i)=tanh(s)
        pe(i)=-tanh(s*bl/pl)+1
        bze(i)=sqrt(2*(1+.5*(bze0**2-bye(i)**2)-pe(i)))
        jze(i)=(1-(tanh(s))**2)/bl
        vxs(i)=0.
        ! vys(i)=vy0*cos(0.25*pi*s*bl)
        ! vys(i)=vy0*tanh(s*bl/lv)
        vys(i)=vy0
        ! vys(i)=vy0*(-1./(cosh(s*bl/lv))+1)
        vzs(i)=0.
        rhoe(i)=1.
        do j=1,nj
            x(i,j,6)=psie(i)
            x(i,j,7)=bze(i)
            x(i,j,1)=rhoe(i)
            x(i,j,2)=pe(i)
            x(i,j,3)=vxs(i)
            x(i,j,4)=vys(i)
            x(i,j,5)=vzs(i)
        enddo
    enddo
    

    do i=1,ni
        s=(i-(ni+1)/2)*dx/bl
        sii=exp(-s*s)*am
        do j=1,nj
            ! sij=sin(mm*2.*pi*(j-1)*dy/4.)
            sij=-cos(mm*pi*(j-1)*dy)+1
            ! sij=sin(0.75*(j-1)*dy*pi)
            ! sij=sin(mm*2.*pi*(j-1)*dy/4.)
            seed(i,j)=sij*sii
            x(i,j,6)=x(i,j,6)+seed(i,j)
            ! x(i,j,6)=0.+sij*sii
        enddo
    enddo
end Subroutine initial

!********************************************************************************
