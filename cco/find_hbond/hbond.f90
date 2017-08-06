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

      bPUMP(2) = bPUMP(1)
      bCHEM(2) = bCHEM(1)
      bCheckedWater = .false.
      if((nGLU.eq.3) .and. (.not. (bPUMP(2) .and. bCHEM(2)))) call check_all_hbond_from_glu(bPUMP(2),bCHEM(2))

      bPUMP(3) = bPUMP(2)
      bCHEM(3) = bCHEM(2)
      bCheckedWater = .false.
      if(.not. (bPUMP(3) .and. bCHEM(3))) call check_all_hbond(bPUMP(3),bCHEM(3))
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
!      do k=1,nneigh
!        idx = k*3 - 2 + nRef
!        do l=1,nGLU
!!          distMat(l,idx,:) = tempDistMat(idxGLU(l),idxWatNeigh(k),:)
!          distMat(l,idx,:) = tempDistMat(l,idx,:)
!          distMat(idx,l,:) = distMat(l,idx,:)
!        enddo
!        do l=nGLU+1,nGLU+nPRD
!!          distMat(l,idx,:) = tempDistMat(idxPRD(l-nGLU),idxWatNeigh(k),:)
!          distMat(l,idx,:) = tempDistMat(l,idx,:)
!          distMat(idx,l,:) = distMat(l,idx,:)
!        enddo
!!        do l=nGLU+nPRD+1,nGLU+nPRD+nBNC
!!          distMat(l,k,:) = distMat(traj % pos(:,idxBNC(l-nGLU-nPRD)),traj % pos(:,idxWatNeigh(k)),:)
!!          distMat(k,l,:) = distMat(l,k,:)
!!        enddo
!      enddo
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
        atidx = idxWatNeigh(temp) + mod((idx-nRef-1),3) 
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
          distMat(j,i,:3) = -distMat(i,j,:3)
          distMat(j,i,4) = distMat(i,j,4)
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
        if(dOH(4).gt.dOO(4)) return
        call get_dist(idxH,idxO2,dHO)
        if(dHO(4).gt.dOO(4)) return
        angle = dot_product(dHO(:3),dOH(:3))/(dHO(4)*dOH(4))
        angle = dacos(angle)
!        bHbond = .true.
        if(angle.le.angHbond) then
          bHbond = .true.
!          write(*,*) angle,angHbond,'angle'
!          write(*,*) dOO(:3),dOH(:3),dHO(:3)
!          write(*,*) dOO(4),dOH(4),dHO(4)
        endif
      endif

    end subroutine check_hbond

    ! subroutine check_cub(idxO,bLigand) : check whether water is ligated to CuB
    subroutine check_cub(idxO,bLigand)
      use variables
      implicit none
      integer, intent(in) :: idxO
      logical, intent(out) :: bLigand
      integer :: idxCu
      real*8,dimension(4) :: dOCu

      idxCu = nGLU + nPRD + 1
      bLigand = .false.

      call get_dist(idxO,idxCu,dOCu)
      if(dOCu(4).le.dCuB) bLigand = .true.

    end subroutine check_cub

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
      real*8,dimension(4) :: dOO, dOH, dHO
      integer :: atidx1,atidx2

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
!        write(*,*) 'neigh search',idxO,idx
        call check_hbond(idxO,idxH,idx,bHbond)
        if(bHbond) then
          call check_acceptor(idx,idx+1,bPUMP,bCHEM)
          if(bPUMP.and.bCHEM) return

          call check_acceptor(idx,idx+2,bPUMP,bCHEM)
          if(bPUMP.and.bCHEM) return
        endif

      enddo

    end subroutine check_acceptor

    ! recursive subroutine check_all(idxO,bPUMP,bCHEM) : check hbond network
    recursive subroutine check_all(idxO,bPUMP,bCHEM)
      use variables
      implicit none
      integer, intent(in) :: idxO
      logical, intent(inout) :: bPUMP, bCHEM

      integer :: i,j,idx,cnt,idxWat            ! counters
      integer :: idxH1, idxH2
      logical :: bHbond
      real*8              :: angle
      real*8,dimension(3) :: vOH, vOO
      real*8,dimension(4) :: dOO, dOH

      if(bPUMP.and.bCHEM) return

      idxWat = (idxO+2-nRef)/3
      if(idxWat.gt.0) bCheckedWater(idxWat) = .true.

      ! check for H-bond to PRD
      i=nGLU+1
      call check_hbond(idxO,idxO+1,i,bHbond)
      if(.not. bHbond) call check_hbond(idxO,idxO+2,i,bHbond)
      if(.not. bHbond) call check_hbond(idxO,idxO+1,i+1,bHbond)
      if(.not. bHbond) call check_hbond(idxO,idxO+2,i+1,bHbond)
      if((.not. bHbond) .and. (nPRD .eq. 3)) call check_hbond(idxO,i+2,i+1,bHbond)
      if(bHbond) bPUMP = .true.
      if(bPUMP.and.bCHEM) return

      ! check for H-bond to BNC-OH or BNC-H2O
      i=nGLU+nPRD+2
      call check_hbond(idxO,idxO+1,i,bHbond)
      if(.not. bHbond) call check_hbond(idxO,idxO+2,i,bHbond)
      if(.not. bHbond) call check_hbond(idxO,i+1,i,bHbond)
      if((.not. bHbond) .and. (nBNC.eq.4)) call check_hbond(idxO,i+2,i,bHbond)
      ! check for H-bond to CuB
      if((.not. bHbond)) call check_cub(idxO,bHbond)
      if(bHbond) bCHEM = .true.
      if(bPUMP.and.bCHEM) return

      do i=1,nneigh
        if(bCheckedWater(i)) cycle  ! pass if already checked

        idx = i*3 - 2 + nRef
        call check_hbond(idxO,idxO+1,idx,bHbond)
        if(.not. bHbond) call check_hbond(idxO,idxO+2,idx,bHbond)
        if(.not. bHbond) call check_hbond(idxO,idx+1,idx,bHbond)
        if(.not. bHbond) call check_hbond(idxO,idx+2,idx,bHbond)
        if(bHbond) then
          call check_all(idx,bPUMP,bCHEM)
          if(bPUMP.and.bCHEM) return
        endif

      enddo

    end subroutine check_all

    ! recursive subroutine check_all_hbond(bPUMP,bCHEM) : check hbond network
    recursive subroutine check_all_hbond(bPUMP,bCHEM)
      use variables
      implicit none
      logical, intent(inout) :: bPUMP, bCHEM

      integer :: i,idx,idxWat            ! counters
      logical :: bHbond

      if(bPUMP.and.bCHEM) return

      do i=1,nneigh
        if(bCheckedWater(i)) cycle  ! pass if already checked

        idx = i*3 - 2 + nRef
        call check_hbond(1,idx+1,idx,bHbond)
        if(.not. bHbond)  call check_hbond(1,idx+2,idx,bHbond)
        if(.not. bHbond)  call check_hbond(2,idx+1,idx,bHbond)
        if(.not. bHbond)  call check_hbond(2,idx+2,idx,bHbond)
        if((.not. bHbond) .and. (nGLU.eq.3)) call check_hbond(2,3,idx,bHbond)
        if(bHbond) then
          call check_all(idx,bPUMP,bCHEM)
          if(bPUMP.and.bCHEM) return
        endif

      enddo

    end subroutine check_all_hbond

    ! subroutine check_all_hbond_from_glu(bPUMP,bCHEM) : check hbond network
    subroutine check_all_hbond_from_glu(bPUMP,bCHEM)
      use variables
      implicit none
      logical, intent(inout) :: bPUMP, bCHEM

      integer :: i,idx,idxWat            ! counters
      logical :: bHbond

      if(bPUMP.and.bCHEM) return

      if(nGLU.ne.3) return
      do i=1,nneigh
        if(bCheckedWater(i)) cycle  ! pass if already checked

        idx = i*3 - 2 + nRef
        call check_hbond(2,3,idx,bHbond)
        if(bHbond) then
          call check_all(idx,bPUMP,bCHEM)
          if(bPUMP.and.bCHEM) return
        endif

      enddo

    end subroutine check_all_hbond_from_glu

end module hbond
