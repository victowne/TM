Program mhd2d
    USE omp_lib
    use variable
    implicit none
    integer cause,it,i,j,ttt
    real*8 t
    
    open(15,file='BXM.TXT',status='replace')
    open(14,file='EQUILIBRIUM.TXT',status='replace')
    open(13,file='INPSI.TXT',status='replace')
    open(12,file='BPSIM.TXT',status='replace')
    open(11,file='BWIDTH.TXT',status='replace')
    open(10,file='ENERGY.TXT',status='replace')
    open(9,file='JMIN.TXT',status='replace')
    open(8,file='BVELOCITY.TXT',status='replace')
    open(7,file='BFORCE.TXT',status='replace')
    open(70,file='BRATE.TXT',status='replace')

    call initial
    do i=1,ni
        write(14,'(6ES18.8)') psie(i), bye(i), pe(i), bze(i), rhoe(i), jze(i)
    end do
    write(13,'(801ES18.8)') ((seed(i,j),i=1,ni),j=1,nj)

    !main loop
    ttt=0
    do it=1,nt
        ! write(*,*)it
        tt=dt*it

        if(mod(it,nplot1)==0) then
        ttt=ttt+1
            call plot1(ttt)
            ! call plot2(ttt)
        endif
        
        if(mod(it-1,nplot2)==0) then
            ! call bxmax
            ! write(15,'(F10.4,ES18.8)') tt, log(bxm)

            call psimaxwidth
            write(11,'(F10.4,ES18.8)') tt, width*dx
            write(8,'(F10.4,2I6)') tt, psixmax, psiymax
            write(70,'(F10.4,ES18.8)') tt, psimin

            ! call calcenergy(x)
            ! write(10,'(F10.4,4ES18.8)') tt, em, ek, eh, et

            call calcj(x)
            write(12,'(F10.4,ES18.8)') tt, psimin
            write(9,'(F10.4,ES18.8)') tt, eta*jmin

            call calcf(x)
            write(7,'(F10.4,4ES18.8)') tt,fp,fmx,fmy,fmz
        end if

        t=t+dt
        call rk4

        do i=1,ni
            do j=1,nj
                cause=1
                if(x(i,j,1)<0) exit
                cause=2
                if(x(i,j,2)<0) exit
            enddo
        enddo
        cause=3
        if(kw==1)       exit
        cause=0

    enddo
    !main loop end
    call exitinfo(cause,it)

    close(15)
end Program mhd2d

!********************************************************************************
Subroutine calcbxy(xi)                          ! calculate Bx, Bz, need check!!!
    USE omp_lib
    use variable
    implicit none
    integer i,j,k,jp1,jm1
    real*8, dimension(ni,nj,7) :: xi     ! 1 2 3 4 5 6 7 --> rho p ux uy uz psi bz
    
    ! calculate Bx, By, need check!!!
    CALL omp_set_num_threads(numthreads)
    !$omp parallel do private(i,j) 
    do i=2,ni-1
        do j=1,nj
            jp1=j+1
            if(j==nj) jp1=1
            jm1=j-1
            if(j==1) jm1=nj
            bx(i,j)=dy2*(xi(i,jp1,6)-xi(i,jm1,6))       ! rdx=0.5/dx, rdz=0.5/dz
            by(i,j)=-dx2*(xi(i+1,j,6)-xi(i-1,j,6))
        enddo
    enddo
    !$omp end parallel do 

    ! by boundary
    do j=1,nj
        jp1=j+1
        if(j==nj) jp1=1
        jm1=j-1
        if(j==1) jm1=nj
        ! bx(1,j)=-rdz*(xi(1,jp1,5)-xi(1,jm1,5))
        by(1,j)=dx2*(3.*xi(1,J,6)-4.*xi(2,J,6)+xi(3,J,6))
        by(ni,j)=dx2*(-xi(ni-2,J,6)+4.*xi(ni-1,J,6)-3.*xi(ni,J,6))
    enddo
    ! bx boundary
    do j=1,nj
        jp1=j+1
        if(j==nj) jp1=1
        jm1=j-1
        if(j==1) jm1=nj
        bx(1,j)=dy2*(xi(1,jp1,6)-xi(1,jm1,6))
        bx(ni,j)=dy2*(xi(ni,jp1,6)-xi(ni,jm1,6))
    enddo

end Subroutine calcbxy

!********************************************************************************
Subroutine bxmax                                               ! calculate Bx max

    use variable
    implicit none
    integer i,j
    real*8 sss
    
    call calcbxy(x)                                                ! Calculate Bx
    
    bxm=-100.0
    do i=1,ni
        do j=1,nj
            sss=abs(bx(i,j))
            if(sss>bxm) bxm=sss
        enddo
    enddo
    if(bxm<small) bxm=small

end Subroutine bxmax
!********************************************************************************
Subroutine calcj(xi)                          ! calculate Bx, Bz, need check!!!
    USE omp_lib
    use variable
    implicit none
    integer i,j,k,jp1,jm1
    real*8, dimension(ni,nj,7) :: xi     ! 1 2 3 4 5 6 7 --> rho p ux uy uz psi bz
    
    real*8, dimension(ni,nj) :: ppsi
    real*8, dimension(ni) :: ppsimin
    ! real*8  num(1)
    ! integer  xmax,ymin
    real*8  num(2)
    integer  jxmax,jymax
    ! psi(:,:)=xi(:,:,6)

    ! psimax=-99999999
    ! psimin=99999999
    ! width=0.
    ! do i=1,ni
    !       do j=1,nj
    !               if(psimax<psi(i,j))then
    !               psimax=psi(i,j)
    !               psixmax=i
    !               psiymax=j
    !               end if
    !       enddo
    ! enddo

    ! ppsimin(:)=psi(:,(nj+1)/2)
    ! psimin=maxval(ppsimin)
    ! num=maxloc(ppsimin)
    ! xmax=num(1)
    ! calculate Jx, Jy, Jz
    do i=1,ni
        do j=1,nj
            jp1=j+1
            if(j==nj) jp1=1
            jm1=j-1
            if(j==1) jm1=nj
            jx(i,j)=dy2*(xi(i,jp1,7)-xi(i,jm1,7))
        enddo
    enddo
    CALL omp_set_num_threads(numthreads)
    !$omp parallel do private(i,j) 
    do i=2,ni-1
        do j=1,nj
            jp1=j+1
            if(j==nj) jp1=1
            jm1=j-1
            if(j==1) jm1=nj
            jy(i,j)=-dx2*(xi(i+1,j,7)-xi(i-1,j,7))
            jz(i,j)=-ddx*(xi(i+1,j,6)+xi(i-1,j,6)-2.*xi(i,j,6))&
                    -ddy**(xi(i,jp1,6)+xi(i,jm1,6)-2.*xi(i,j,6))
            jy(1,j)=jy(2,j)
            jz(1,j)=jz(2,j)
            jy(ni,j)=jy(ni-1,j)
            jz(ni,j)=jz(ni-1,j)
        enddo
    enddo
   !$omp end parallel do 

    ! jmin=jz(xmax,(nj+1)/2)
    jmin=maxval(jz)
    num=maxloc(jz)
    jxmax=num(1)
    jymax=num(2)
    psimin=xi(jxmax,jymax,6)
end Subroutine calcj
                        
!********************************************************************************
Subroutine calcjb(xi)                          
    USE omp_lib
    use variable
    implicit none
    integer i,j,k,jp1,jm1,ip1,im1
    real*8, dimension(ni,nj,7) :: xi     ! 1 2 3 4 5 6 7 --> rho p ux uy uz psi bz

    CALL omp_set_num_threads(numthreads)
    !$omp parallel do private(i,j) 
    do j=1,nj
        do i=2,ni-1
            im1=i-1
            ip1=i+1
            jb(i,j)=jbs*dx2*(xi(ip1,j,2)-xi(im1,j,2))
        enddo
        jb(1,j)=jb(2,j)
        jb(ni,j)=jb(ni-1,j)
    enddo
   !$omp end parallel do 

end Subroutine calcjb
!********************************************************************************
Subroutine calcf(xi)                          
    USE omp_lib
    use variable
    implicit none
    integer i,j,k,jp1,jm1,ip1,im1
    real*8, dimension(ni,nj,7) :: xi     ! 1 2 3 4 5 6 7 --> rho p ux uy uz psi bz
    real*8  numm(2)
    integer  psixmin,psiymin

    CALL omp_set_num_threads(numthreads)
    !$omp parallel do private(i,j) 
    do j=1,nj
        do i=1,ni
            psi(i,j)=xi(i,j,6)
        end do 
    end do 
    !$omp end parallel do 
    numm=minloc(psi)
    psixmin=numm(1)
    psiymin=numm(2)
    fp=dx2*(xi(psixmin+1,psiymin,2)-xi(psixmin-1,psiymin,2))&
    & +dy2*(xi(psixmin,psiymin+1,2)-xi(psixmin,psiymin-1,2))
    
    call calcj(xi)
    call calcbxy(xi)
    fmx=jy(psixmin,psiymin)*xi(psixmin,psiymin,7)-jz(psixmin,psiymin)*by(psixmin,psiymin)
    fmy=jz(psixmin,psiymin)*bx(psixmin,psiymin)-jx(psixmin,psiymin)*xi(psixmin,psiymin,7)
    fmz=jx(psixmin,psiymin)*by(psixmin,psiymin)-jy(psixmin,psiymin)*bx(psixmin,psiymin)

end Subroutine calcf

!********************************************************************************
Subroutine psimaxwidth                                               ! calculate psimax

    use variable
    implicit none
    integer i,j

    call rk4                                               
    psi(:,:)=x(:,:,6)

    psimax=-99999999
    psimin=99999999
    width=0.
    do i=1,ni
        do j=1,nj
            if(psimax<psi(i,j))then
                psimax=psi(i,j)
                psixmax=i
                psiymax=j
            end if
        enddo
    enddo

    do j=1,nj
        psimin=min(psimin,psi(psixmax,j))
    end do

    do i=2,ni
        if(psi(i,psiymax)>psimin)then
            if(width<0.000001)then
                width=(psi(i,psiymax)-psimin)/(psi(i,psiymax)-psi(i-1,psiymax))
            else
                width=width+1
            end if
        else if(width>0.000001)then
            width=width+(psi(i-1,psiymax)-psimin)/(psi(i-1,psiymax)-psi(i,psiymax))
            exit
        end if
    end do
end Subroutine psimaxwidth
!********************************************************************************
Subroutine calcenergy(xi)                          ! calculate Bx, Bz, need check!!!
    USE omp_lib
    use variable
    implicit none
    integer i,j,k,jp1,jm1
    real*8, dimension(ni,nj,7) :: xi     ! 1 2 3 4 5 6 7 --> rho p ux uy uz psi bz
    
    call calcbxy(xi)

    em=0.
    ek=0.
    eh=0.
    CALL omp_set_num_threads(numthreads)
    !$omp parallel do private(i,j) 
    do j=1,nj
        do i=1,ni
            em=em+.5*(bx(i,j)**2+by(i,j)**2+xi(i,j,7)**2)
            ek=ek+.5*xi(i,j,1)*(xi(i,j,3)**2+xi(i,j,4)**2+xi(i,j,5)**2)
            eh=eh+xi(i,j,2)/(gamma-1)
            et=em+ek+eh
        enddo
    enddo
    !$omp end parallel do
end Subroutine calcenergy

!********************************************************************************
