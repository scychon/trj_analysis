!****************************************************************
! 
! this module finds the h-bonded molecules within cutoff distance
! The search is done for the molecules in the array
!
!****************************************************************

module hbond


contains
 
    ! perform all hbond calculation
    subroutine find_hbond(bPUMP,bCHEM)
      use variables
      implicit none
      integer :: i,j,k,idx              ! counters
      logical,dimension(2), intent(out) :: bPUMP, bCHEM

      allocate(bCheckedWater(nneigh))
      bPUMP = .false.
      bCHEM = .false.
      bCheckedWater = .false.

!      write(*,*) 'Shrinking the distMat'
      call shrink_distMat
!      write(*,*) 'Calculating the distMat'
      call calc_distMat

!      write(*,*) 'Checking acceptor network'
      if(nGLU.eq.3) call check_acceptor(2,3,bPUMP(1),bCHEM(1))
!      call find_donor(bPUMP(1),bCHEM(1))
      deallocate(bCheckedWater)
      
    end subroutine find_hbond

    ! subroutine shrink_distMat : Shrink the distance matrix to save memory access time
    subroutine shrink_distMat
      use variables
      implicit none
      integer :: k,l,idx              ! counters

      if(allocated(distMat)) deallocate(distMat)
      allocate(distMat(nAtCavity,nAtCavity,4))
      distMat = 0
      do k=1,nneigh
        idx = k*3 - 2 + nRef
        do l=1,nGLU
!          distMat(l,idx,:) = tempDistMat(idxGLU(l),idxWatNeigh(k),:)
          distMat(l,idx,:) = tempDistMat(l,idx,:)
          distMat(idx,l,:) = distMat(l,idx,:)
        enddo
        do l=nGLU+1,nGLU+nPRD
!          distMat(l,idx,:) = tempDistMat(idxPRD(l-nGLU),idxWatNeigh(k),:)
          distMat(l,idx,:) = tempDistMat(l,idx,:)
          distMat(idx,l,:) = distMat(l,idx,:)
        enddo
!        do l=nGLU+nPRD+1,nGLU+nPRD+nBNC
!          distMat(l,k,:) = distMat(traj % pos(:,idxBNC(l-nGLU-nPRD)),traj % pos(:,idxWatNeigh(k)),:)
!          distMat(k,l,:) = distMat(l,k,:)
!        enddo
      enddo
!      deallocate(tempDistMat)
       
    end subroutine shrink_distMat

    ! subroutine get_dist(idx1,idx2,dist) : get the distance between two different idxs in distMat
    subroutine get_dist(idx1,idx2,dist)
      use variables
      implicit none
      integer, intent(in) :: idx1,idx2              ! counters
      real*8,dimension(4),intent(out) :: dist
      integer :: atidx1,atidx2              ! counters

      dist(:) = distMat(idx1,idx2,:)
      if(dist(4).eq.0) then
        call get_atIdx(idx1,atidx1)
        call get_atIdx(idx2,atidx2)
        dist(:3) = traj % pos(:,atidx2) - traj % pos(:,atidx1)
        dist(4) = sqrt(dot_product(dist(:3),dist(:3)))
      endif

    end subroutine get_dist

    ! subroutine get_atIdx(idx, atidx) : get the atom idx from idx in distMat
    subroutine get_atIdx(idx, atidx)
      use variables
      implicit none
      integer, intent(in) :: idx              ! idx in water neighborlist
      integer, intent(out) :: atidx           ! idx in atom index
      integer :: temp                         ! counters

      if(idx.le.nGLU) then
        atidx = idxGLU(idx)
      elseif(idx.le.(nGLU+nPRD)) then
        atidx = idxPRD(idx-nGLU)
      elseif(idx.le.nRef) then
        atidx = idxBNC(idx-nPRD-nGLU)
      else
        temp = (idx-nRef-1)/3 + 1
        atidx = idxWatNeigh(temp) + mod((idx-nRef),3) -1
      endif

    end subroutine get_atIdx

    ! subroutine calc_distMat : calculate complete distance matrix for
    ! neighborlist
    subroutine calc_distMat
      use variables
      implicit none
      integer :: i,j              ! counters

!      write(*,*) 'nAtCavity : ', nAtCavity
      do i=1,nAtCavity
        do j=i+1,nAtCavity
          call get_dist(i,j,distMat(i,j,:))
          distMat(j,i,:) = distMat(i,j,:)
        enddo
      enddo

    end subroutine calc_distMat

    ! subroutine check_hbond(idxO,idxH,idxO2,bHbond) : check whether OH-O combi is hbond or not
    subroutine check_hbond(idxO,idxH,idxO2,bHbond)
      use variables
      implicit none
      integer, intent(in) :: idxO, idxH, idxO2
      logical, intent(out) :: bHbond
      real*8              :: angle
      real*8,dimension(4) :: dOO, dOH, dHO

      bHbond = .false.

      call get_dist(idxO,idxO2,dOO)
      if(dOO(4).le.dHbond) then
        call get_dist(idxO,idxH,dOH)
        call get_dist(idxH,idxO2,dHO)
        angle = dot_product(dHO(:3),dOH(:3))/(dHO(4)*dOH(4))
        angle = dacos(angle)
        if(angle.le.angHbond) bHbond = .true.
      endif

    end subroutine check_hbond

    ! recursive subroutine check_acceptor(idxO,idxH,bPUMP,bCHEM) : check hbond network
    recursive subroutine check_acceptor(idxO,idxH,bPUMP,bCHEM)
      use variables
      implicit none
      integer, intent(in) :: idxO, idxH
      logical, intent(inout) :: bPUMP, bCHEM

      integer :: i,j,idx,cnt,idxWat            ! counters
      logical :: bHbond
      real*8              :: angle
      real*8,dimension(3) :: vOH, vOO
      real*8,dimension(4) :: dOO, dOH

      if(bPUMP.and.bCHEM) return

      idxWat = (idxO+2-nRef)/3
      if(idxWat.gt.0) bCheckedWater(idxWat) = .true.

      do i=nGLU+1,nGLU+2
        call check_hbond(idxO,idxH,i,bHbond)
        if(bHbond) bPUMP = .true.
      enddo
      if(bPUMP.and.bCHEM) return

      do i=nGLU+nPRD+1,nGLU+nPRD+2
        call check_hbond(idxO,idxH,i,bHbond)
        if(bHbond) bCHEM = .true.
      enddo
      if(bPUMP.and.bCHEM) return

      do i=1,nneigh
        if(bCheckedWater(i)) cycle  ! pass if already checked

        idx = i*3 - 2 + nRef
        call check_hbond(idxO,idxH,idx,bHbond)
        if(bHbond) then
          call check_acceptor(idx,idx+1,bPUMP,bCHEM)
          if(bPUMP.and.bCHEM) return

          call check_acceptor(idx,idx+2,bPUMP,bCHEM)
          if(bPUMP.and.bCHEM) return
        endif

      enddo

    end subroutine check_acceptor

end module hbond
