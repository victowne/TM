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
