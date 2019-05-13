            
Module variable

	implicit none
	
	! grid parameters
	integer, parameter :: ni=801                      ! nx, x grid number
	integer, parameter :: nj=201                      ! ny, y grid number
	integer, parameter :: nt=25000000
	real*8, parameter :: dx=0.005, dy=0.02
	real*8 :: dt=0.002    ! giving a dt < min(dx,dz)/[sqrt(1.0+0.5*gamma*beta)*va]
	real*8 :: dx2=0.5/dx, dy2=0.5/dy, ddx=1.0/(dx*dx), ddy=1.0/(dy*dy)
	real*8 :: tt 
	! physical parameters
	real*8, parameter :: gamma=1.66667              ! parameters in mhd equations
	real*8, parameter :: eta=0.00005
	real*8, parameter :: vis=0.0001
	real*8, parameter :: dp=0.0001
	real*8, parameter :: drho=0.0001
	real*8, parameter :: beta=0.5                                  ! beta in x=Lx
	real*8, parameter :: rho0=1.0             ! rho0 -- mass density rho0 in x=Lx
	real*8, parameter :: b0=1.0                                        ! b0 -- B0
	real*8, parameter :: bze0=3.0                                      ! bze0 
	real*8, parameter :: bl=0.8                    ! bl -- width of current sheet
	real*8, parameter :: pl=0.6                 ! bl -- width of initial pressure
	real*8 :: va                                                ! Alfven velocity
	real*8 :: p0                                      ! p0 -- pressure p0 in x=Lx
	real*8 :: t0                                 ! t0 -- temperature T0 in x=Lx
	real*8 :: bxm                                           ! bxm -- log(max(Bx))
	real*8 :: em, ek, eh, et
	
	! other parameters
	integer, parameter :: nplot1=50000, nplot2=500, nplot3=500000
	real*8, parameter :: pi=3.14159265358979, small=1.0e-10, large=1.0e10
	real*8, parameter :: cfl=1.5! Courant�CFriedrichs�CLewy condition parameter
	real*8, parameter :: am=0.0001, mm=0.5                    ! am -- prtb amplitude
	real*8, parameter :: abd=1.0                ! abd -- perturbation width in x
	integer ::  kw=0                     ! kw -- warning, =0 normal, =1 divergent

	real*8, dimension(ni,nj,7) :: x    ! 1 2 3 4 5 6 7 --> rho p ux uy uy psi bz
	real*8, dimension(ni,nj,7) :: y, f1, f2, f3, f4         ! temp arrays for rk4
	real*8, dimension(ni,nj) :: bx, by, pt, jb                           ! pt=p+b^2/2
	real*8, dimension(ni,nj) :: jx, jy, jz
	real*8, dimension(ni) :: psie, bye, pe, bze, rhoe, jze                !equilibrium

	real*8, dimension(ni,nj) ::psi                          !for psimax and width
	real*8 :: psimax,psimin,width
	real*8 :: fp, fmx, fmy, fmz
	integer :: psixmax,psiymax 
	real*8 jmin
	real*8, dimension(ni,nj) ::vor                         !vorticiy
	real*8, dimension(ni,nj) ::seed
	real*8, dimension(ni,nj) ::fvx, fvy, fvz, fjbx, fjby, fjbz, fh
    
	!forced cancel
	real*8 :: ta, tau, tt0, ttrmp
	real*8, parameter :: numta=200, numta1=1000, numta2=800
	real*8, parameter :: psip=0.1, pper=0., pp0=0.02, f0=0.

	!neoclassical current
	real*8, parameter :: jbs=0.08

	!shear flow
	real*8, dimension(ni) :: vxs, vys, vzs                                   
	real*8, parameter :: vx0=0.0, vy0=0.07, vz0=0.0      
	real*8, parameter :: lv=0.8

	!rmp
	real*8, parameter :: widthon=9000                                  !start rmp
	real*8, parameter :: psiper=0.3                        
	real*8 :: phi=1.0*pi                                    !the phase diffrence


	!eccd
	real*8, parameter :: eccdon=19500 
	real*8, parameter :: ejz=0.0, delta=0.6
	real*8 :: delta2=delta**2
 
	integer, parameter :: numthreads=24
	
end Module variable

!********************************************************************************
Program mhd2d
	USE omp_lib
	use variable
	implicit none
	integer cause,it,i,j,ttt
	INTEGER(4) :: time_begin, time_end, time_rate
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
	call system_clock(time_begin,time_rate) 
	ttt=0
	do it=1,nt
	! write(*,*)it
	tt=dt*it

		if(mod(it,nplot1)==0) then
		ttt=ttt+1
			call plot1(ttt)
			! call plot2(ttt)
		endif
		
		if(mod(it,nplot3)==0) then
		ttt=ttt+1
			call plot2(ttt)
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
				if(x(i,j,1)<0)	exit
				cause=2
				if(x(i,j,2)<0)	exit
			enddo
		enddo
		cause=3
		if(kw==1)	exit
		cause=0
		
	enddo
	!main loop end
	call exitinfo(cause,it)

	close(15)
	CALL system_clock(time_end,time_rate)
	WRITE(*,*) 'time is: ',(time_end - time_begin)/time_rate	
end Program mhd2d

!********************************************************************************
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
Subroutine rk4	                              ! 4-th Runge-Kutta time integration
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
		! 	xi(ni,j,6)=psiper*sin(mm*2.*(j-1)*dy*pi/4+phi)*tanh((tt-widthon)/ttrmp)!rmp
		! 	! xi(ni,j,6)=-psip*(1-exp(-tt/tau))*sin(mm*2.*pi*(j-1)*dy/4.)
		! else 
		! 	! xi(ni,j,6)=-psip*(1-exp(-tt/tau))*sin(mm*2.*pi*(j-1)*dy/4.)!force
		! 	xi(ni,j,6)=2.*xi(ni-1,j,6)-xi(ni-2,j,6)
		! end if  
		! xi(ni,j,6)=psip*tanh(tt/tau)*sin(mm*2.*pi*(j-1)*dy/4.)*(1-tanh((tt-tt0)/tau))
		! xi(ni,j,6)=psip*tanh(tt/tau)*sin(mm*2.*pi*(j-1)*dy/4.)
		xi(ni,j,7)=xi(ni-1,j,7)
	end do
	return
	end Subroutine boundary

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
	! 	do j=1,nj
	! 		if(psimax<psi(i,j))then
	! 		psimax=psi(i,j)
	! 		psixmax=i
	! 		psiymax=j
	! 		end if
	! 	enddo
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
Subroutine calcftotal(xi)                          
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

	call calcj(xi)
	call calcbxy(xi)
	!$omp parallel do private(i,j) 
	do j=1,nj
		jp1=j+1
		jm1=j-1
		if(j==nj) jp1=1
		if(j==1) jm1=nj
		do i=2,ni-1
		ip1=i+1
		im1=i-1
			fvx(i,j)=vis*(ddx*(xi(ip1,j,3)+xi(im1,j,3)-2.0*xi(i,j,3))+ddy*(xi(i,jp1,3)+xi(i,jm1,3)-2.0*xi(i,j,3)))
			fvy(i,j)=vis*(ddx*(xi(ip1,j,4)+xi(im1,j,4)-2.0*xi(i,j,4))+ddy*(xi(i,jp1,4)+xi(i,jm1,4)-2.0*xi(i,j,4)))
			fvz(i,j)=vis*(ddx*(xi(ip1,j,5)+xi(im1,j,5)-2.0*xi(i,j,5))+ddy*(xi(i,jp1,5)+xi(i,jm1,5)-2.0*xi(i,j,5)))
			fjbx(i,j)=jy(i,j)*xi(i,j,7)-jz(i,j)*by(i,j)
			fjby(i,j)=jz(i,j)*bx(i,j)-jx(i,j)*xi(i,j,7)
			fjbz(i,j)=jx(i,j)*by(i,j)-jy(i,j)*bx(i,j)
			fh(i,j)=dx2*(xi(ip1,j,2)-xi(im1,j,2))+dy2*(xi(i,jp1,2)-xi(i,jm1,2))
		end do 
	end do 
	!$omp end parallel do 
	do j=1,nj
		fvx(1,j)=fvx(2,j)
		fvy(1,j)=fvy(2,j)
		fvz(1,j)=fvz(2,j)
		fjbx(1,j)=fjbx(2,j)
		fjby(1,j)=fjby(2,j)
		fjbz(1,j)=fjbz(2,j)
		fh(1,j)=fh(2,j)

		fvx(ni,j)=fvx(ni-1,j)
		fvy(ni,j)=fvy(ni-1,j)
		fvz(ni,j)=fvz(ni-1,j)
		fjbx(ni,j)=fjbx(ni-1,j)
		fjby(ni,j)=fjby(ni-1,j)
		fjbz(ni,j)=fjbz(ni-1,j)
		fh(ni,j)=fh(ni-1,j)
	end do
end Subroutine calcftotal

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
Subroutine calcvorcity(xi)                          ! calculate Bx, Bz, need check!!!
	USE omp_lib
	use variable
	implicit none
	integer i,j,k,jp1,jm1
	real*8, dimension(ni,nj,7) :: xi     ! 1 2 3 4 5 6 7 --> rho p ux uy uz psi bz
	CALL omp_set_num_threads(numthreads)
	!$omp parallel do private(i,j)
	do j=1,nj
		jp1=j+1
		if(j==nj) jp1=1
		jm1=j-1
		if(j==1) jm1=nj
			do i=2,ni-1
			vor(i,j)=dx2*(xi(i+1,j,4)-xi(i-1,j,4)) &
				-dy2*(xi(i,jp1,3)-xi(i,jm1,3))
			enddo
			vor(1,j)=vor(2,j)
			vor(ni,j)=vor(ni-1,j)
	enddo
	!$end omp parallel do private(i,j)
end Subroutine calcvorcity
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
SUBROUTINE smooth(xi)
	USE omp_lib	
	use variable
	implicit none
	integer i,j,k,jp1,jm1
	real*8, dimension(ni,nj,7) :: xi 
	real*8, dimension(ni,nj) :: w
	real*8  cf
	cf=0.99
	CALL omp_set_num_threads(numthreads)
    do K=1,5
		!$omp parallel do private(i,j)
		do j=1,nj
		jp1=j+1
		if(j==nj) jp1=1
		jm1=j-1
		if(j==1) jm1=nj
			do i=2,ni-1
			W(I,J)=xi(I,J,K)+(1.-CF)*0.25*( &
				(xi(I-1,J,K)-2.*xi(I,J,K)+xi(I+1,J,K))&
				+(xi(I,JM1,K)-2.*xi(I,J,K)+xi(I,JP1,K)))
			end do
			w(ni,j)=xi(ni,j,k)+(1.-cf)*0.25*(xi(ni,j,k)+xi(ni-1,j,k)&
			+xi(ni,jm1,k)+xi(ni,jp1,k))
		end do        
        !$end omp parallel do 

		!$omp parallel do private(i,j)
		DO  J=1,NJ
			DO  I=2,NI-1
			xi(I,J,K)=W(I,J)
			end do
		end do
		!$omp end parallel do
    END DO

	k=7
		!$omp parallel do private(i,j)
		do j=1,nj
		jp1=j+1
		if(j==nj) jp1=1
		jm1=j-1
		if(j==1) jm1=nj
			do i=2,ni-1
			W(I,J)=xi(I,J,K)+(1.-CF)*0.25*( &
				(xi(I-1,J,K)-2.*xi(I,J,K)+xi(I+1,J,K))&
				+(xi(I,JM1,K)-2.*xi(I,J,K)+xi(I,JP1,K)))
			end do
			w(ni,j)=xi(ni,j,k)+(1.-cf)*0.25*(xi(ni,j,k)+xi(ni-1,j,k)&
			+xi(ni,jm1,k)+xi(ni,jp1,k))
		end do   
		!$omp end parallel do 
          
		DO  J=1,NJ
			DO  I=2,NI-1
			xi(I,J,K)=W(I,J)
			end do
		end do
	RETURN 
END SUBROUTINE smooth

!********************************************************************************
Subroutine exitinfo(cause,it)

	use variable
	implicit none

	integer cause, it
	select case(cause)
		case(1)
			write(*,*) 'Negative density at',it,'-th step'
		case(2)
			write(*,*) 'Negative perssure at',it,'-th step'
		case(3)
			write(*,*) 'Divergent at',it,'-th step'
		case default
			write(*,*) 'it = nt_end=',it,'-th step'
	end select
end Subroutine exitinfo

!********************************************************************************
subroutine plot1(ttt)         ! Routine used to write contour data to external file

	use variable
	character(len=80) :: filename1,filename2,filename3
	character(len=80) :: filename4,filename5,filename6
	character(len=80) :: filename7,filename8,filename9
	character(len=80) :: filename10,filename11,filename12,filename13,filename14
	character(len=80) :: form1,form2,form3,form4,form5,form6,form7,form8,form9
	character(len=80) :: form10,form11,form12,form13,form14

	integer i,j,ttt

	write(form1,*)ttt
	write(form2,*)ttt
	write(form3,*)ttt
	write(form4,*)ttt
	write(form5,*)ttt
	write(form6,*)ttt
	write(form7,*)ttt
	write(form8,*)ttt
	write(form9,*)ttt
	write(form10,*)ttt
	write(form11,*)ttt
	write(form12,*)ttt
	! write(form13,*)ttt
	write(form14,*)ttt

	filename1="NTMRHO"//trim(adjustl(form1))//".TXT"
	filename2="NTMP"//trim(adjustl(form2))//".TXT"
	filename3="NTMVX"//trim(adjustl(form3))//".TXT"
	filename4="NTMVY"//trim(adjustl(form4))//".TXT"
	filename5="NTMVZ"//trim(adjustl(form5))//".TXT"
	filename6="NTMPSI"//trim(adjustl(form6))//".TXT"
	filename7="NTMBZ"//trim(adjustl(form7))//".TXT"
	! filename8="NTMBX"//trim(adjustl(form8))//".TXT"
	! filename9="NTMBY"//trim(adjustl(form9))//".TXT"
	! filename10="NTMJX"//trim(adjustl(form10))//".TXT"
	! filename11="NTMJY"//trim(adjustl(form11))//".TXT"
	filename12="NTMJZ"//trim(adjustl(form12))//".TXT"
	! filename13="NTMVOR"//trim(adjustl(form13))//".TXT"
	filename14="NTMJB"//trim(adjustl(form14))//".TXT"

	OPEN(16,file=filename1)
	OPEN(17,file=filename2)
	OPEN(18,file=filename3)
	OPEN(19,file=filename4)
	OPEN(20,file=filename5)
	OPEN(21,file=filename6)
	OPEN(22,file=filename7)
	! OPEN(23,file=filename8)
	! OPEN(24,file=filename9)
	! OPEN(25,file=filename10)
	! OPEN(26,file=filename11)
	OPEN(27,file=filename12)
	! OPEN(28,file=filename13)
	OPEN(29,file=filename14)

	call calcbxy(x)   
	call calcj(x) 
	call calcvorcity(x) 
	call calcjb(x)                                     
	
	write(16,'(801E15.6)') (((x(i,j,1)),i=1,ni),j=1,nj)
	write(17,'(801E15.6)') (((x(i,j,2)),i=1,ni),j=1,nj)
    write(18,'(801E15.6)') (((x(i,j,3)),i=1,ni),j=1,nj)
    write(19,'(801E15.6)') (((x(i,j,4)),i=1,ni),j=1,nj)
    write(20,'(801E15.6)') (((x(i,j,5)),i=1,ni),j=1,nj)
	write(21,'(801E15.6)') (((x(i,j,6)),i=1,ni),j=1,nj)
	write(22,'(801E15.6)') (((x(i,j,7)),i=1,ni),j=1,nj)
	! write(23,'(801E15.6)') (((bx(i,j)),i=1,ni),j=1,nj)
	! write(24,'(801E15.6)') (((by(i,j)),i=1,ni),j=1,nj)
	! write(25,'(801E15.6)') (((jx(i,j)),i=1,ni),j=1,nj)
	! write(26,'(801E15.6)') (((jy(i,j)),i=1,ni),j=1,nj)
	write(27,'(801E15.6)') (((jz(i,j)),i=1,ni),j=1,nj)
	! write(28,'(801E15.6)') (((vor(i,j)),i=1,ni),j=1,nj)
	write(29,'(801E15.6)') (((jb(i,j)),i=1,ni),j=1,nj)
end subroutine plot1
!********************************************************************************
subroutine plot2(ttt)         ! Routine used to write contour data to external file

	use variable
	character(len=80) :: filename15,filename16,filename17,filename18,filename19,filename20,filename21
	character(len=80) :: form15,form16,form17,form18,form19,form20,form21
	integer i,j,ttt

	write(form15,*)ttt
	write(form16,*)ttt
	write(form17,*)ttt
	write(form18,*)ttt
	write(form19,*)ttt
	write(form20,*)ttt
	write(form21,*)ttt

	! filename15="NTMFVX"//trim(adjustl(form15))//".TXT"
	! filename16="NTMFVY"//trim(adjustl(form16))//".TXT"
	filename17="NTMFVZ"//trim(adjustl(form17))//".TXT"
	! filename18="NTMFJBX"//trim(adjustl(form18))//".TXT"
	! filename19="NTMFJBY"//trim(adjustl(form19))//".TXT"
	filename20="NTMFJBZ"//trim(adjustl(form20))//".TXT"
	filename21="NTMFH"//trim(adjustl(form21))//".TXT"

	! OPEN(30,file=filename15)
	! OPEN(31,file=filename16)
	OPEN(32,file=filename17)
	! OPEN(33,file=filename18)
	! OPEN(34,file=filename19)
	OPEN(35,file=filename20)
	OPEN(36,file=filename21)

	call calcftotal(x)                                        
	
	! write(30,'(801E15.6)') (((fvx(i,j)),i=1,ni),j=1,nj)
	! write(31,'(801E15.6)') (((fvy(i,j)),i=1,ni),j=1,nj)
	write(32,'(801E15.6)') (((fvz(i,j)),i=1,ni),j=1,nj)
	! write(33,'(801E15.6)') (((fjbx(i,j)),i=1,ni),j=1,nj)
	! write(34,'(801E15.6)') (((fjby(i,j)),i=1,ni),j=1,nj)	
	write(35,'(801E15.6)') (((fjbz(i,j)),i=1,ni),j=1,nj)	
	write(36,'(801E15.6)') (((fh(i,j)),i=1,ni),j=1,nj)
end subroutine plot2
