!****************************************************************
! 
! this module reads topology file and save it to the type topfile
!
!****************************************************************

module calc_corr


contains
 
    ! calculate correlation function
    subroutine calc_corr_IP
      use variables
      implicit none
      integer :: i,j,k,idx              ! counters
      logical :: tf

      allocate(corr_IP(nxtc))
      allocate(corr_IP_cat(nxtc))
      allocate(corr_IP_ani(nxtc))

      write(*,*) 'start calculating corr_IP'
      OPEN (UNIT=6,FORM='FORMATTED')
      corr_IP = 0
      corr_IP_cat = 0
      corr_IP_ani = 0
      do i=1,nxtc
        do j=i,nxtc
          do idx=1,nmolcat
            tf = idxPairIon(j,idx) .eq. idxPairIon(i,idx)
            if(tf) corr_IP_cat(j-i+1) = corr_IP_cat(j-i+1) + 1
            if((i==j) .and. (.not. tf)) write(*,*) i,j,idx,idxPairIon(j,idx),idxPairIon(i,idx)
          enddo
          do idx=nmolcat+1,nmolcat+nmolani
            tf = idxPairIon(j,idx) .eq. idxPairIon(i,idx)
            if(tf) corr_IP_ani(j-i+1) = corr_IP_ani(j-i+1) + 1
            if((i==j) .and. (.not. tf)) write(*,*) i,j,idx,idxPairIon(j,idx),idxPairIon(i,idx)
          enddo
        enddo
        corr_IP = corr_IP_cat + corr_IP_ani
        write(6,100,advance='no') achar(13),i,'th frame has finished  '
      enddo
      do i=1,nxtc
        corr_IP_cat(i) = corr_IP_cat(i)/(nmolcat*(nxtc-i+1))
        corr_IP_ani(i) = corr_IP_ani(i)/(nmolani*(nxtc-i+1))
        corr_IP(i) = corr_IP(i)/(nmolsys*(nxtc-i+1))
      enddo
      write(*,*) 'calculating corr_IP finished'
 100    FORMAT(A,I5,A)
      close(6)

    end subroutine calc_corr_IP


    subroutine calc_corr_CP
      use variables
      implicit none
      integer :: i,j,k,idx              ! counters
      logical :: tf

      allocate(corr_CP(nxtc))
      allocate(corr_CP_cat(nxtc))
      allocate(corr_CP_ani(nxtc))
      write(*,*) 'start calculating corr_CP'

      OPEN (UNIT=6,FORM='FORMATTED')
      corr_CP = 0
      corr_CP_cat = 0
      corr_CP_ani = 0
      write(*,*) nmolcat, 'cathions and ',nmolani,' anions'
      do i=1,nxtc
        do j=i,nxtc
          do idx=1,nmolcat
            tf = ALL(idxClustIons(j,idx,:) .eq. idxClustIons(i,idx,:))
            if(tf) corr_CP_cat(j-i+1) = corr_CP_cat(j-i+1) + 1
          enddo
          do idx=nmolcat+1,nmolcat+nmolani
            tf = ALL(idxClustIons(j,idx,:) .eq. idxClustIons(i,idx,:))
            if(tf) corr_CP_ani(j-i+1) = corr_CP_ani(j-i+1) + 1
          enddo
        enddo
        write(6,100,advance='no') achar(13),i,'th frame has finished  '
      enddo
      corr_CP = corr_CP_cat + corr_CP_ani

      do i=1,nxtc
        corr_CP_cat(i) = corr_CP_cat(i)/(nmolcat*(nxtc-i+1))
        corr_CP_ani(i) = corr_CP_ani(i)/(nmolani*(nxtc-i+1))
        corr_CP(i) = corr_CP(i)/(nmolsys*(nxtc-i+1))
      enddo
      write(*,*) 'calculating corr_CP finished'
 100    FORMAT(A,I5,A)
      close(6)

    end subroutine calc_corr_CP

end module calc_corr
