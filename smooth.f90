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
END SUBROUTINE smooth

!********************************************************************************
