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
