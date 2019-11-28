Module variable

    implicit none
    
    ! grid parameters
    integer, parameter :: ni=801                      ! nx, x grid number
    integer, parameter :: nj=201                      ! ny, y grid number
    integer, parameter :: nt=2500000                  ! steps
    real*8, parameter :: dx=0.005, dy=0.02            ! deltax, deltay
    real*8 :: dt=0.002                                ! giving a dt < min(dx,dz)/[sqrt(1.0+0.5*gamma*beta)*va]
    real*8 :: dx2=0.5/dx, dy2=0.5/dy, ddx=1.0/(dx*dx), ddy=1.0/(dy*dy)
    real*8 :: tt 
    ! physical parameters
    real*8, parameter :: gamma=1.66667                
    real*8, parameter :: eta=0.00005                  ! resisitivity
    real*8, parameter :: vis=0.0001                   ! viscosity
    real*8, parameter :: dp=0.0001              
    real*8, parameter :: drho=0.0001
    real*8, parameter :: beta=0.5                     ! beta in x=Lx
    real*8, parameter :: rho0=1.0                     ! rho0 -- mass density rho0 in x=Lx
    real*8, parameter :: b0=1.0                       ! b0 -- B0
    real*8, parameter :: bze0=3.0                     ! bze0 
    real*8, parameter :: bl=0.2                       ! bl -- width of current sheet
    real*8, parameter :: pl=0.6                       ! bl -- width of initial pressure
    real*8 :: va                                      ! Alfven velocity
    real*8 :: p0                                      ! p0 -- pressure p0 in x=Lx
    real*8 :: t0                                      ! t0 -- temperature T0 in x=Lx
    real*8 :: bxm                                     ! bxm -- log(max(Bx))
    real*8 :: em, ek, eh, et
    
    ! other parameters
    ! output interval
    integer, parameter :: nplot1=25000, nplot2=10000 
    real*8, parameter :: pi=3.14159265358979, small=1.0e-10, large=1.0e10
    real*8, parameter :: cfl=1.5                      ! CFL condition parameter
    real*8, parameter :: am=0.0001, mm=0.5            ! am -- prtb amplitude
    real*8, parameter :: abd=1.0                      ! abd -- perturbation width in x
    integer ::  kw=0                                  ! kw -- warning, =0 normal, =1 divergent

    real*8, dimension(ni,nj,7) :: x                   ! 1 2 3 4 5 6 7 --> rho p ux uy uy psi bz
    real*8, dimension(ni,nj,7) :: y, f1, f2, f3, f4   ! temp arrays for rk4
    real*8, dimension(ni,nj) :: bx, by, pt, jb        ! pt=p+b^2/2
    real*8, dimension(ni,nj) :: jx, jy, jz
    real*8, dimension(ni) :: psie, bye, pe, bze, rhoe, jze                !equilibrium

    real*8, dimension(ni,nj) ::psi                    ! for psimax and width
    real*8 :: psimax,psimin,width
    real*8 :: fp, fmx, fmy, fmz
    integer :: psixmax,psiymax 
    real*8 jmin
    real*8, dimension(ni,nj) ::vor                    ! vorticiy
    real*8, dimension(ni,nj) ::seed
    real*8, dimension(ni,nj) ::fvx, fvy, fvz, fjbx, fjby, fjbz, fh

    !forced cancel
    real*8 :: ta, tau, tt0, ttrmp
    real*8, parameter :: numta=200, numta1=1000, numta2=800
    real*8, parameter :: psip=0.1, pper=0., pp0=0.02, f0=0.

    !neoclassical current
    real*8, parameter :: jbs=0.00

    !shear flow
    real*8, dimension(ni) :: vxs, vys, vzs                                   
    real*8, parameter :: vx0=0.0, vy0=0.01, vz0=0.0      
    real*8, parameter :: lv=0.8

    !rmp
    real*8, parameter :: widthon=9000                  ! start rmp
    real*8, parameter :: psiper=0.3                        
    real*8 :: phi=1.0*pi                               ! the phase diffrence

    !eccd
    real*8, parameter :: eccdon=1.d8 
    real*8, parameter :: ejz=0.0, delta=0.6
    real*8 :: delta2=delta**2

    integer, parameter :: numthreads=24
        
end Module variable

!********************************************************************************
