!  XDR Fortran Interface xtc Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!
 
program calc_corr_time
 
    ! 1. Use the xdr interface
    use xtc, only: xtcfile
    use topol, only: topfile
    use variables
    use calc_corr
 
    implicit none

    ! 2. Declare a variable of type xtcfile
    type(xtcfile) :: traj
    type(topfile) :: sys

!  --------------------------------------------------------------------------
  integer            :: nbins=200                   !# of bins 
  integer            :: nsize=1000                  !# of bins 
  integer            :: i,j,k,l,idx,ixtc,idxmol,idxatom           !counters
  integer            :: isize,ifile                  !counters
  integer            :: narg, cptArg, nxtcfile       !#of arg & counter of arg


  logical::lookForInp=.FALSE.
  logical::lookForOut=.FALSE.
  logical::fileExist

  real*8             :: delr=0.02d0,rmax=3.0d0,boxsize=0.3d0      !
  real*8             :: mass, molmass, dr, rclustcat, rclustani
  real*8,dimension(3)    :: xpos,sysbox,compos,dist,box
  real*8,dimension(300)  :: rbin
  real*8, allocatable,dimension(:,:,:) :: comtraj,comtrajcat,comtrajani,temp   ! trajectory of com of each molecules (frame idx, mol idx, xyz)
  real*8, allocatable,dimension(:,:) :: boxtraj,tempbox      !matrix for # of water molecules along z direction
  real*8, allocatable,dimension(:) :: temptime     !matrix for # of water molecules along z direction
  character(len=256)         :: strout,strtitle
  character(len=256)         :: strInFile, strOutFile, strTopFile
  character(len=256)         :: strListFile, strIPFile, strCPFile
  character(len=256)         :: str, line
  character(len=256),allocatable,dimension(:) :: strXtcfiles
! ----------------------------------------------------------------


  rclustcat = 0.73d0
  rclustani = 0.73d0
! Initialization
  strInFile = ""
  strOutFile = ""
  strTopFile = 'param_bmimbf4.dat'
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
    read(7, *) strTopFile
    read(7, '(A)') line
    read(7, *) strListFile
    read(7, *) strIPFile
    read(7, *) strCPFile
    close(7)

    ! 3. Read the topology file and initialize atomic mass data for system
    write(*,'(A,A)') 'reading topoly file ',trim(strTopFile)
    write(*,*) ''
    call sys % init(strTopFile)



    ! 3. Initialize it with the name of traj file you want to read in.
    ixtc=0
    do ifile=1,nxtcfile
    strInFile = trim(strXtcfiles(ifile))
    write(*,'(A,A)') 'reading trajectory file ',trim(strInFile)
    write(*,*) ''
    call traj % init(strInFile)

    nmolsys = sys % numsysmol
    allocate(comtraj(nsize,nmolsys,3))
    nmolcat = sys % moltype(1) % nummol
    nmolani = sys % moltype(2) % nummol

    allocate(timestamp(nsize))
    allocate(boxtraj(nsize,3))
 
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
        ixtc = ixtc + 1

        ! record current time
        timestamp(ixtc) = traj % time

        ! check the size of nwat matrix
        isize = size(comtraj(:,1,1))
        if(isize .lt. ixtc) then
          allocate(temp(isize+nsize,nmolsys,3))
          temp=0
          temp(:isize,:,:)=comtraj
          deallocate(comtraj)
          allocate(comtraj(isize+nsize,nmolsys,3))
          comtraj=temp
          deallocate(temp)

          allocate(tempbox(isize+nsize,3))
          tempbox=0
          tempbox(:isize,:)=boxtraj
          deallocate(boxtraj)
          allocate(boxtraj(isize+nsize,3))
          boxtraj=tempbox
          deallocate(tempbox)

          allocate(temptime(isize+nsize))
          temptime=0
          temptime(:isize)=timestamp
          deallocate(timestamp)
          allocate(timestamp(isize+nsize))
          timestamp=temptime
          deallocate(temptime)

          nsize = nsize*2
        endif

        do k=1,3
          boxtraj(ixtc,k) = traj % box(k,k)
        enddo

        idxatom = 0
        idxmol = 0
        do i=1, sys % nummoltype
          do j=1, sys % moltype(i) % nummol
            idxmol = idxmol + 1
            compos = 0

            do k=1, sys % moltype(i) % numatom
              idxatom = idxatom +1
              mass = sys % moltype(i) % atommass(k)
              xpos(:) = traj % pos(:,idxatom)
              compos(:) = compos(:) + xpos(:)*mass
            enddo

            molmass = sys % moltype(i) % molarmass
            compos(:) = compos(:)/molmass
            comtraj(ixtc,idxmol,:) = compos(:)
!            write(*,*) comtraj(ixtc,idxmol,:)
          enddo
        enddo 

        write(6,100) ixtc,'th frame has finished  ' 
 100    FORMAT('+', I5,A)  

        call traj % read
 
    end do

    ! 5. Close the file
    call traj % close
    end do

    nxtc = ixtc

    write(6,*) 'start analyzing the data' 

    ! check the size of trajectory matrix
    isize = size(comtraj(:,1,1))
    if(isize .gt. nxtc) then
      allocate(temp(nxtc,nmolsys,3))
      temp=0
      temp(:,:,:)=comtraj(:nxtc,:,:)
      deallocate(comtraj)
      allocate(comtraj(nxtc,nmolsys,3))
      comtraj=temp
      deallocate(temp)
    endif

    write(6,*) 'copying position data into each ion matrices' 
    
    allocate(comtrajcat(nxtc,nmolcat,3))
    allocate(comtrajani(nxtc,nmolani,3))
    allocate(idxPairIon(nxtc,nmolsys))
    allocate(distPairIon(nxtc,nmolsys))
    allocate(idxClustIons(nxtc,nmolsys,10))
    comtrajcat(:,:,:) = comtraj(:,:nmolcat,:)
    comtrajani(:,:,:) = comtraj(:,(nmolcat+1):,:)

    write(6,*) 'start generating neighborlist data'
    idx = 0
    write(6,*) 'first frame start'
    do ixtc=1, nxtc
      box(:) = boxtraj(ixtc,:)
      dr = norm(box*0.5d0)
      distPairIon(ixtc,:) = dr
      do i=1,nmolcat
        do j=1,nmolani
          dist(:) = comtrajani(ixtc,j,:) - comtrajcat(ixtc,i,:)
          dr = getdr(dist, box)
          if(dr .lt. 0.1) then
            write(*,*) ixtc, i, j, dist,comtrajani(ixtc,j,:),comtrajcat(ixtc,i,:)
            write(*,'(6f15.6)') comtraj(ixtc,i,:),comtrajcat(ixtc,i,:)
            write(*,'(6f15.6)') comtraj(ixtc,nmolcat+j,:),comtrajani(ixtc,j,:)
          endif

          ! Check if anion is closest to cathion
          if(dr .lt. rclustcat) then
            do k=1, 10
              if(idxClustIons(ixtc,i,k) .eq. 0) then
                idxClustIons(ixtc,i,k) = nmolcat+j
                exit
              endif
            enddo
            if(distPairIon(ixtc,i) .gt. dr) then
              idxPairIon(ixtc,i) = nmolcat+j
              distPairIon(ixtc,i) = dr
            endif
          endif

          ! Check if cathion is closest to anion
          if(dr .lt. rclustani) then
            do k=1, 10
              if(idxClustIons(ixtc,nmolcat+j,k) .eq. 0) then
                idxClustIons(ixtc,nmolcat+j,k) = i
                exit
              endif
            enddo
            if(distPairIon(ixtc,nmolcat+j) .gt. dr) then
              idxPairIon(ixtc,nmolcat+j) = i
              distPairIon(ixtc,nmolcat+j) = dr
            endif
          endif

        enddo 
      enddo 
      write(6,100) ixtc,'th frame has finished  ' 
    enddo

    close(6)

    open(18,file=strListFile)
    do i=1,nxtc
      write(18,'(11I10)') i, idxPairIon(i,1:nmolsys:40)
    enddo
    close(18) 


    open(16,file=strIPFile)
    call calc_corr_IP
    do i=1,nxtc
      write(16,'(I5,3f15.9)') i, corr_IP(i), corr_IP_cat(i), corr_IP_ani(i)
    enddo
    close(16) 

    open(17,file=strCPFile)
    call calc_corr_CP
    do i=1,nxtc
      write(17,'(I5,3f15.9)') i, corr_CP(i), corr_CP_cat(i), corr_CP_ani(i)
    enddo
    close(17) 


contains

real*8 function getdr( vec_dr, vec_box )
    implicit none
    real*8, intent(in), dimension(3) :: vec_box, vec_dr
    real*8, dimension(3) :: tempvec

    tempvec(:) = vec_dr(:) - (nint(vec_dr(:)/vec_box(:)))*vec_box(:)
    getdr = norm(tempvec)
end function getdr

real*8 function norm( vec )
    implicit none
    real*8, intent(in), dimension(3) :: vec

    norm = sqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))

end function norm

end program calc_corr_time
