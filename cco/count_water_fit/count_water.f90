!  XDR Fortran Interface xtc Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!
 
program read_xtc_prog
 
    ! 1. Use the xdr interface
    use xtc, only: xtcfile
 
    implicit none

    ! 2. Declare a variable of type xtcfile
    type(xtcfile) :: traj
!  --------------------------------------------------------------------------
  integer            :: nbins=200                   !# of bins 
  integer            :: nsize=1000                  !# of bins 
  integer            :: i,j,k,l,idx,ixtc,ifile      !counters
  integer            :: idxatom                     !counters
  integer            :: nxtc,isize                  !counters
  integer            :: narg, cptArg, nxtcfile      !#of arg & counter of arg
  integer            :: nCentRef      !#of reference coordinates for the proten center
  integer            :: idxwatb, idxwate    ! index of first water and last water
  integer            :: nframe     ! average per every nframe trajectories

  integer, dimension(3) :: ibox,nbox
  integer, allocatable,dimension(:) :: isInternalWater   ! 1 if the atom is an internal water oxygen, 0 otherewise
  integer, allocatable,dimension(:) :: idxcenteratoms   !idxs of reference atoms for protein center
  integer, allocatable,dimension(:,:,:,:) :: nwat,temp   !matrix for # of water molecules
  integer, allocatable,dimension(:,:) :: nwat_z,temp_z      !matrix for # of water molecules along z direction
  real, allocatable,dimension(:,:) :: nwatcore,nwatd,nwatupcore       !matrix for statics of # of water molecules in key sections(core, d-channel, upper area of core
  real, allocatable,dimension(:,:,:,:) :: nwatstat       !matrix for statics of # of water molecules
! nwatstat : 1:average, 2:stdv, 3:first value, 4:last value, 5:min 6:max
  real, allocatable,dimension(:) :: zrefs   !idxs of reference atoms for protein center
  real, allocatable,dimension(:,:) :: xcavityrefs   !positions of reference atoms for internal cavity

  logical::lookForInp=.FALSE.
  logical::lookForOut=.FALSE.
  logical::fileExist

  real*8             :: delr=0.02d0,rmax=3.0d0,boxsize=0.3d0      !
  real*8             :: nwatavg
  real*8             :: rCavitySize
  real*8,dimension(3)    :: xpos,xcenter
  real*8,dimension(3)    :: sysbox,boxunit,sizeprot
  real*8,dimension(3)    :: xshift          ! relertive position of reference atoms in protein
  real*8,dimension(300)  :: rbin
!  real*8, allocatable,dimension(:) :: re_to_e_list,temp  ! e_to_e data list 
  real*8             :: r, re_to_e, rsg_N, rsg_N_eig1,rsg_N_eig2,rsg_N_eig3
  character(len=256)         :: strout,strtitle
  character(len=256)         :: strInFile, strOutFile,strZmatFile,strTotmatFile
  character(len=256)         :: str, line
  character(len=256),allocatable,dimension(:) :: strXtcfiles
! ----------------------------------------------------------------

  boxsize = 0.25d0
  sysbox = (/ 11.2, 11.2, 13.6 /)
  nbox(:) = sysbox(:)/boxsize +1
!  nbox = (/ 40, 40, 47 /)

! Initialization
  strInFile = ""
  strOutFile = ""
!Check if any arguments are found
 narg=command_argument_count()
!Loop over the arguments
 if(narg>0)then
!loop across options
  do cptArg=1,narg
   call get_command_argument(cptArg,str)
   select case(adjustl(str))
     case("--help","-h")
      write(*,*)"This is program TestArg : Version 0.1"

!First known args
    case("-f")
     lookForInp=.TRUE. !change logical value
    case("-o")
     lookForOut=.TRUE.

    case default
!Treat the second arg of a serie
     if(LookForInp)then
      strInFile=adjustl(str) !assign a value to pedfile
      LookForInp=.FALSE. !put the logical variable to its initial value
     elseif(LookForOut)then
      strOutFile=adjustl(str)
      inquire(file=strOutFile,exist=fileExist)
      if(fileExist)then
       write(*,*)'file ',strOutFile,' exist'
      endif
      LookForOut=.FALSE.
     else
      write(*,*)"Option ",adjustl(str),"unknown"
     endif
    end select
  end do
 end if

if(strInFile .eq. "") then
  write(*,*) 'input file has not set'
  write(*,*) 'usage : count_water -f infile.xtc -o outfile.dat -n index.ndx'
  stop
endif

inquire(file=strInFile,exist=fileExist)!check if it exist
if(.not.fileExist)then
  write(*,*)'file ',strInFile,' not found'
  write(*,*) 'usage : count_water -f infile.xtc -o outfile.dat -n index.ndx'
  stop
endif

if(strOutFile .eq. "") then
  write(*,*) 'usage : count_water -f infile.xtc -o outfile.dat -n index.ndx'
  write(*,*) 'output file has not set'
  write(*,*) 'will use default outfile name matWater.dat'
  strOutFile = 'matWater.dat'
endif

    OPEN(unit=10,file='zref.dat')
    ! Read param file to locate strXtcfiles
    OPEN(unit=7,file=strInFile,status='old')
    read(7, '(A)') line
    write(*,*) trim(line)
    read(7,*) nxtcfile
    write(*,*) nxtcfile
    allocate(strXtcfiles(nxtcfile))
    do ifile=1,nxtcfile
      read(7,'(A)') strXtcfiles(ifile)
      write(*,*) trim(strXtcfiles(ifile))
    enddo
    read(7, '(A)') line
    read(7, *) boxunit
    read(7, '(A)') line
    read(7, *) sizeprot
    read(7, '(A)') line
    read(7, *) xshift
    read(7, '(A)') line
    read(7, *) idxwatb, idxwate
    read(7, '(A)') line
    read(7, *) nCentRef
    allocate(idxcenteratoms(nCentRef))
    allocate(zrefs(nCentRef))
    allocate(xcavityrefs(3,nCentRef))
    read(7, *) idxcenteratoms
    read(7, '(A)') line
    read(7, *) rCavitySize
    read(7, '(A)') line
    read(7, *) strOutFile
    read(7, *) strZmatFile
    read(7, *) strTotmatFile
    read(7, '(A)') line
    read(7, *) nframe
    close(7)

    nbox(:) = sizeprot(:)/boxunit(:) +1

! read in the input parameters
!    open(15,file=strInFile)
!    do i=1,5
!      read(15,'(A)') strtitle    ! skip header
!      write(*,*)strtitle
!      read(15,'(A)') strout    ! skip header
!      write(*,*)strout
!      if(trim(strout).eq.'@ s4 legend "<R\sg\N> eig3"') exit
!    enddo    

!   ----------------------------------------------------------------
!   write title

!          if(size(r_list) .lt. j) then
!            allocate(temp(size(r_list)+nsize))
!            temp=0
!            temp(:size(r_list))=r_list
!            deallocate(r_list)
!            allocate(r_list(size(temp)))
!            r_list=temp
!            deallocate(temp)
!
!            allocate(temp_E(3,size(r_list)))
!            temp_E=0
!            temp_E(:,:size(E_list,2))=E_list(:,:)
!            deallocate(E_list)
!            allocate(E_list(3,size(temp_E,2)))
!            E_list=temp_E
!            deallocate(temp_E)
!
!          endif

    ! 3. Initialize it with the name of traj file you want to read in.
    allocate(nwat(nbox(1),nbox(2),nbox(3),nsize))
    nwat=0
    allocate(nwat_z(nbox(3),nsize))
    nwat_z=0
    allocate(nwatstat(nbox(1),nbox(2),nbox(3),6))
    nwatstat=0

    ixtc=0
    do ifile=1,nxtcfile
    strInFile = trim(strXtcfiles(ifile))
    write(*,'(A,A)') 'reading trajectory file ',trim(strInFile)
    write(*,*) ''
    call traj % init(strInFile)

!    nmolsys = sys % numsysmol
!    allocate(comtraj(nsize,nmolsys,3))
!    nmolcat = sys % moltype(1) % nummol
!    nmolani = sys % moltype(2) % nummol

!    allocate(timestamp(nsize))
!    allocate(boxtraj(nsize,3))

 
    ! 4. Read in each configuration. Everything is stored in the xtcfile
    ! type (precision, time,
    !    step, no of atoms, positions, etc.). Look in the xtc module for
    !    more details.
    !    You can save the positions in the loop for your calculations in
    !    another array, or 
    !    do your calculations after each read.
 
    call traj % read

    OPEN (UNIT=6,FORM='FORMATTED',CARRIAGECONTROL='FORTRAN')
        ! Just an example to show what was read in
        write(6,'(a,f12.6,a,i0)') " Time (ps): ", traj % time, "  Step: ", traj % STEP
        write(6,'(a,f12.6,a,i0)') " Precision: ", traj % prec, "  No. Atoms: ", traj % NATOMS
        ! This is the same order as found in the GRO format fyi
        write(6,'(9f9.5)')  traj % box(1,1), traj % box(2,2), traj % box(3,3), &
                            traj % box(1,2), traj % box(1,3), & 
                            traj % box(2,1), traj % box(2,3), &
                            traj % box(3,1), traj % box(3,2)
        write(6,*) '' 
    do while ( traj % STAT == 0 )

        ixtc= ixtc+1

        ! check the size of nwat matrix
        isize = size(nwat(1,1,1,:))
        if(isize .lt. ixtc) then
          allocate(temp(nbox(1),nbox(2),nbox(3),isize+nsize))
          temp=0
          temp(:,:,:,:isize)=nwat
          deallocate(nwat)
          allocate(nwat(nbox(1),nbox(2),nbox(3),isize+nsize))
          nwat=temp
          deallocate(temp)

          allocate(temp_z(nbox(3),isize+nsize))
          temp_z=0
          temp_z(:,:isize)=nwat_z
          deallocate(nwat_z)
          allocate(nwat_z(nbox(3),isize+nsize))
          nwat_z=temp_z
          deallocate(temp_z)

          nsize = nsize*2
        endif

        xcenter = 0
        do j=1,nCentRef
           idxatom = idxcenteratoms(j)
           xcenter(:) = xcenter(:) + traj % pos(:,idxatom)
           xcavityrefs(:,j) = xcavityrefs(:,j) + traj      % pos(:,idxatom)
           zrefs(j) = traj % pos(3,idxatom)
        enddo
        xcenter(:) = xcenter(:)/nCentRef - xshift(:)
        write(10,'(4f12.6)') zrefs,xcenter(3)
 
!        write(*,'(3f9.3)') traj % pos
        do j=idxwatb,idxwate,3
          xpos(:) = traj % pos(:,j) - xcenter(:)
          ibox(:) = xpos(:)/boxunit(:)+1
!          write(*,'(3f9.3)') xpos
          if( (ibox(1).gt.0) .and. (ibox(1).lt.nbox(1)) .and. (ibox(2).gt.0) .and. &
              (ibox(2).lt.nbox(2)) .and. (ibox(3).gt.0) .and. (ibox(3).lt.nbox(3))) then
            nwat(ibox(1),ibox(2),ibox(3),ixtc) = nwat(ibox(1),ibox(2),ibox(3),ixtc)+1
            nwat_z(ibox(3),ixtc) = nwat_z(ibox(3),ixtc)+1
!            write(*,'(3f9.3,4I3)') xpos, nwat(ibox(1),ibox(2),ibox(3),ixtc),ibox
          endif
          
        enddo


        write(6,100) ixtc,'th frame has finished  ' , nwat_z(4,ixtc)
 100    FORMAT('+', I7,A,I)  

        call traj % read
 
    end do

    ! 5. Close the file
    call traj % close
    enddo

    nxtc = ixtc
    do ixtc=1,nxtc
        nwatstat(:,:,:,1) = nwatstat(:,:,:,1)+nwat(:,:,:,ixtc)
        nwatstat(:,:,:,2) = nwatstat(:,:,:,2)+nwat(:,:,:,ixtc)*nwat(:,:,:,ixtc)
        if (ixtc .eq. 1) then
          nwatstat(:,:,:,3) = nwat(:,:,:,1)
          nwatstat(:,:,:,5) = nwat(:,:,:,1)
          nwatstat(:,:,:,6) = nwat(:,:,:,1)
        else
          do i=1,nbox(1)
            do j=1,nbox(2)
              do k=1,nbox(3)
                if(nwatstat(i,j,k,5) .gt. nwat(i,j,k,ixtc)) nwatstat(i,j,k,5) = nwat(i,j,k,ixtc)
                if(nwatstat(i,j,k,6) .lt. nwat(i,j,k,ixtc)) nwatstat(i,j,k,6) = nwat(i,j,k,ixtc)
              enddo 
            enddo 
          enddo 
          
        endif
    enddo

    nwatstat(:,:,:,4) = nwat(:,:,:,nxtc)
    nwatstat(:,:,:,1) = nwatstat(:,:,:,1)/real(nxtc)
    nwatstat(:,:,:,2) = sqrt(nwatstat(:,:,:,2)/real(nxtc) - nwatstat(:,:,:,1)**2)
    close(10)

    open(16,file=strOutFile)
    open(17,file=strZmatFile)
    open(18,file=strTotmatFile)

    write(str,'(A,I5,A)') '(', nbox(2)*nbox(3),'f9.5)'
    idx = 0
      do i=1,nbox(1)
!        write(*, str) nwatstat(i,:,:,1)
!        write(*, str) nwatstat(i,:,:,2)
        do j=1,nbox(2)
          do k=1,nbox(3)
            if(nwatstat(i,j,k,1) .eq. 0) cycle
!            if((nwatstat(i,j,k,1) .gt. 31.0) .and. (nwatstat(i,j,k,2) .lt. 2.8)) cycle
!            if(nwatstat(i,j,k,6) .le. 2) cycle
            idx = idx + 1
            write(16, '(4I6,7f10.5)') idx, i, j, k, nwatstat(i,j,k,2)/nwatstat(i,j,k,1), nwatstat(i,j,k,:)
          enddo 
        enddo 
      enddo 
    close(16) 

    nwatavg = 0
    do i=1,nbox(1)
      do j=1,nbox(2)
        do k=1,nbox(3)
          idx = 0
          nwatavg = 0
          do ixtc=1,nxtc
             idx = idx + 1
             nwatavg = nwatavg + nwat(i,j,k,ixtc)
             if(mod(idx,nframe) .eq. 0) then
                nwatavg = nwatavg/nframe
                write(18, '(4I6,f10.5)') ixtc, i, j, k, nwatavg
                nwatavg = 0
             endif
          enddo 
        enddo 
      enddo 
    enddo 
    close(18) 

    write(str,'(A,I,A)') '(', nbox(3)+1,'I5)'
    write(6,*) str
    do ixtc=1,nxtc
      write(17, str) ixtc, nwat_z(:,ixtc)
    enddo 
    close(17) 

end program read_xtc_prog
