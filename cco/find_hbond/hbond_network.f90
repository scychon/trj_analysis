!  XDR Fortran Interface xtc Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!
 
program hbond_network
 
    use variables
    use hbond
    
 
    implicit none

!  --------------------------------------------------------------------------
  integer            :: nbins=200                   !# of bins 
  integer            :: nsize=1000                  !# of bins 
  integer            :: nvar =6                     !# of output variables
  integer            :: i,j,k,l,idx,ixtc,ifile      !counters
  integer            :: idxatom                     !counters
  integer            :: isize                       !counters
  integer            :: narg, cptArg, nxtcfile      !#of arg & counter of arg
  integer            :: idxwatb, idxwate    ! index of first water and last water
  integer            :: nframe,iframe     ! average per every nframe trajectories
  integer            :: ineigh     ! counter for water molecules within the cutoff distance
  integer            :: nMaxWaterCnt=200   ! counter for water molecules within the cutoff distance
  integer,dimension(3)     :: nPUMP=0,nCHEM=0                 !resulting total number of PUMP & CHEM conf
  integer,dimension(3)     :: nPUMPframe=0,nCHEMframe=0                 !resulting total number of PUMP & CHEM conf
  integer,allocatable,dimension(:)    :: nframextc
  integer,allocatable,dimension(:,:)    :: nPUMPxtc,nCHEMxtc

  logical::lookForInp=.FALSE.
  logical::lookForOut=.FALSE.
  logical::fileExist
  logical,dimension(2) :: bPUMP, bCHEM

  real*8             :: dist,pi
  real*8,dimension(2)    :: rPUMPxtc,rCHEMxtc
  real*8,dimension(3)    :: xpos, dxglu, dxprd
  real*8,dimension(3)    :: idx_glu,idx_prd
  character(len=256)         :: strout,strtitle
  character(len=256)         :: strInFile, strOutFile,strAvgFile
  character(len=256)         :: str, line
  character(len=256),allocatable,dimension(:) :: strXtcfiles
! ----------------------------------------------------------------


! Initialization
  strInFile = ""
  strOutFile = ""
!Check if any arguments are found
 narg=command_argument_count()
!Loop over the arguments
 pi = dacos(-1.0d0)
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
  write(*,*) 'usage : hbond_network -f param.dat'
  stop
endif

inquire(file=strInFile,exist=fileExist)!check if it exist
if(.not.fileExist)then
  write(*,*)'file ',strInFile,' not found'
  write(*,*) 'usage : hbond_network -f param.dat'
  stop
endif


    ! Read param file to locate strXtcfiles
    OPEN(unit=7,file=strInFile,status='old')
    read(7, '(A)') line
    write(*,*) trim(line)
    read(7,*) nxtcfile
    write(*,*) nxtcfile
    allocate(strXtcfiles(nxtcfile))
    allocate(nPUMPxtc(nxtcfile,3))
    allocate(nCHEMxtc(nxtcfile,3))
    allocate(nframextc(nxtcfile))

    do ifile=1,nxtcfile
      read(7,'(A)') strXtcfiles(ifile)
      write(*,*) trim(strXtcfiles(ifile))
    enddo
    read(7, '(A)') line
    read(7, *) idxGLU
    read(7, '(A)') line
    read(7, *) idxPRD
    read(7, '(A)') line
    read(7, *) idxBNC
    idxCuB = idxBNC(1)
    read(7, '(A)') line
    read(7, *) dCut
    read(7, '(A)') line
    read(7, *) dHbond
!    read(7, '(A)') line
!    read(7, *) dCuB
    dCuB = dHbond
    read(7, '(A)') line
    read(7, *) angHbond
    read(7, '(A)') line
    read(7, *) idxwatb, idxwate
    read(7, '(A)') line
    read(7, *) strOutFile
    read(7, *) strAvgFile
    read(7, '(A)') line
    read(7, *) nframe
    close(7)

    write(*,*) angHbond
!    angHbond = 40
    angHbond = angHbond * pi / 180.0d0
    write(*,*) angHbond
    if(idxGLU(3).gt.0) then
      nGLU = 3
    else
      nGLU = 2
    endif
    if(idxPRD(3).gt.0) then
      nPRD = 3
    else
      nPRD = 2
    endif
    if(idxBNC(4).gt.0) then
      nBNC = 4
    else
      nBNC = 3
    endif
    nRef = nGLU + nPRD + nBNC
    allocate(xglu(nGLU,3))
    allocate(xprd(nPRD,3))
    allocate(xbnc(nBNC,3))

    ! 3. Initialize it with the name of traj file you want to read in.
    allocate(hbondcnt(nvar,nsize))
    hbondcnt=0

    ixtc=0
    do ifile=1,nxtcfile
    strInFile = trim(strXtcfiles(ifile))
    write(*,'(A,A)') 'reading trajectory file ',trim(strInFile)
    write(*,*) ''
    open(16,file=strOutFile)
    open(17,file=strAvgFile)
    write(16,'(A,A)') '# output generated by hbond_network -f ',trim(strInFile)
    write(17,'(A,A)') '# output generated by hbond_network -f ',trim(strInFile)
    call traj % init(strInFile)

 
    ! 4. Read in each configuration. Everything is stored in the xtcfile
    ! type (precision, time,
    !    step, no of atoms, positions, etc.). Look in the xtc module for
    !    more details.
    !    You can save the positions in the loop for your calculations in
    !    another array, or 
    !    do your calculations after each read.
 
    call traj % read
    nPUMPxtc = 0
    nCHEMxtc = 0
    nframextc(ifile) = 0
    

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
        natom = traj % NATOMS
    do while ( traj % STAT == 0 )

        ixtc= ixtc+1
        nframextc(ifile) = nframextc(ifile)+1

        ! check the size of hbondcnt matrix
        isize = size(hbondcnt(1,:))
        if(isize .lt. ixtc) then
          allocate(tempHbond(nvar,isize+nsize))
          tempHbond=0
          tempHbond(:,:isize)=hbondcnt
          deallocate(hbondcnt)
          allocate(hbondcnt(nvar,isize+nsize))
          hbondcnt=tempHbond
          deallocate(tempHbond)

          nsize = nsize*2
        endif

        do j=1,nGLU
          xglu(j,:) = traj % pos(:,idxGLU(j))
!          write(*,*) 'idxGLU', idxGLU(j), xglu(j,:)
        enddo
        do j=1,nPRD
          xprd(j,:) = traj % pos(:,idxPRD(j))
!          write(*,*) 'idxPRD', idxPRD(j), xprd(j,:)
        enddo
        do j=1,nBNC
          xbnc(j,:) = traj % pos(:,idxBNC(j))
!          write(*,*) 'idxBNC', idxBNC(j), xbnc(j,:)
        enddo
        xcub(:) = traj % pos(:,idxCuB)
        xgluc(:) = (xglu(1,:)+xglu(2,:))/2.0d0
        xprdc(:) = (xprd(1,:)+xprd(2,:))/2.0d0
 
!       Find water molecules within the neighborlist cutoff distance
        idx = nRef
        ineigh = 0
!        write(*,*) 'before allocating tempDistMat'

        nneigh = (idxwate-idxwatb)/3+1
        nAtCavity = nneigh*3 + nRef
        
        if(.not. allocated(tempWatNeigh)) allocate(tempWatNeigh(nneigh))
        if(.not. allocated(tempDistMat)) allocate(tempDistMat(nMaxWaterCnt,nMaxWaterCnt,4))
        tempWatNeigh = 0
        tempDistMat = 0
!        write(*,*) 'after allocating tempDistMat'

        ! select only the water molecules within the cutoff distance
        do j=idxwatb,idxwate,3
          idx = idx+1
          xpos(:) = traj % pos(:,j)
          dxglu(:) = abs(xpos(:) - xgluc(:))
          if( (dxglu(1).gt.dCut) .or. (dxglu(2).gt.dCut) .or. (dxglu(3).gt.dCut) ) cycle
!          write(*,*) j,'glu',dxglu(1)
          dxprd(:) = abs(xpos(:) - xprdc(:))
!          write(*,*) 'x',xpos(:), 'glu',xgluc(:),'prd',xprdc(:)
          if( (dxprd(1).gt.dCut) .or. (dxprd(2).gt.dCut) .or. (dxprd(3).gt.dCut) ) cycle

!          write(*,*) j,dxprd(:)
          do k=1,nGLU
            dxglu(:) = xpos(:) - traj % pos(:,idxGLU(k))
            dist = dsqrt(dot_product(dxglu,dxglu))
            if(dist.gt.dCut) goto 150
          enddo
          
          do k=1,nPRD
            dxprd(:) = xpos(:) - traj % pos(:,idxPRD(k))
            dist = dsqrt(dot_product(dxprd,dxprd))
            if(dist.gt.dCut) goto 150
          enddo

          ineigh=ineigh+1
          tempWatNeigh(ineigh) = j
150      enddo

!        write(*,*) 'after water cycle, nneigh :', ineigh
        ! copy temporary idx list into real idxWatNeigh list
        allocate(idxWatNeigh(ineigh))
        idxWatNeigh = tempWatNeigh(:ineigh)
        deallocate(tempWatNeigh)
        nneigh = ineigh

        ! set the number of atoms to consider within the cavity
        nAtCavity = ineigh*3 + nRef

!        write(*,*) 'start finding hbonding network'
        ! calculate hbondcnt network
        call find_hbond(bPUMP,bCHEM)
!        write(*,*) 'found hbonding network'
        deallocate(idxWatNeigh)

        ! hbond acceptor from Glu286H
        if(bPUMP(1)) then
          hbondcnt(1,ixtc)= 1
          nPUMP(1) = nPUMP(1)+1
          nPUMPframe(1) = nPUMPframe(1)+1
          nPUMPxtc(ifile,1) = nPUMPxtc(ifile,1)+1
        ! all hbond link to the acceptor of hbond from Glu286H
        elseif(bPUMP(2)) then
          hbondcnt(3,ixtc)= 1
          nPUMP(2) = nPUMP(2)+1
          nPUMPframe(2) = nPUMPframe(2)+1
          nPUMPxtc(ifile,2) = nPUMPxtc(ifile,2)+1
        ! all hbond link to Glu286
        elseif(bPUMP(3)) then
          hbondcnt(5,ixtc)= 1
          nPUMP(3) = nPUMP(3)+1
          nPUMPframe(3) = nPUMPframe(3)+1
          nPUMPxtc(ifile,3) = nPUMPxtc(ifile,3)+1
        endif

        ! hbond acceptor from Glu286H
        if(bCHEM(1)) then
          hbondcnt(2,ixtc)= 1
          nCHEM(1) = nCHEM(1)+1
          nCHEMframe(1) = nCHEMframe(1)+1
          nCHEMxtc(ifile,1) = nCHEMxtc(ifile,1)+1
        ! all hbond link to the acceptor of hbond from Glu286H
        elseif(bCHEM(2)) then
          hbondcnt(4,ixtc)= 1
          nCHEM(2) = nCHEM(2)+1
          nCHEMframe(2) = nCHEMframe(2)+1
          nCHEMxtc(ifile,2) = nCHEMxtc(ifile,2)+1
        ! all hbond link to Glu286
        elseif(bCHEM(3)) then
          hbondcnt(6,ixtc)= 1
          nCHEM(3) = nCHEM(3)+1
          nCHEMframe(3) = nCHEMframe(3)+1
          nCHEMxtc(ifile,3) = nCHEMxtc(ifile,3)+1
        endif


        write(6,100) ixtc,'th frame has finished  ' , hbondcnt(1,ixtc),hbondcnt(3,ixtc)
 100    FORMAT('+', I7,A,I)  

        call traj % read
        write(16, '(f12.6,I,I,I,I,I,I)') traj % time, hbondcnt(1,ixtc),hbondcnt(1,ixtc)+hbondcnt(3,ixtc),hbondcnt(1,ixtc)+hbondcnt(3,ixtc)+hbondcnt(5,ixtc), hbondcnt(2,ixtc), hbondcnt(2,ixtc)+hbondcnt(4,ixtc),hbondcnt(2,ixtc)+hbondcnt(4,ixtc)+hbondcnt(6,ixtc)
        if(mod(ixtc,nframe).eq.0) then
          nPUMPframe(2) = nPUMPframe(1)+nPUMPframe(2)
          nPUMPframe(3) = nPUMPframe(2)+nPUMPframe(3)
          nCHEMframe(2) = nCHEMframe(1)+nCHEMframe(2)
          nCHEMframe(3) = nCHEMframe(2)+nCHEMframe(3)
          write(17, '(f12.6,I,I,I,f12.6,f12.6,f12.6,I,I,I,f12.6,f12.6,f12.6)') traj % time, nPUMPframe(1),nPUMPframe(2),nPUMPframe(3), real(nPUMPframe(1))/real(nframe),real(nPUMPframe(2))/real(nframe),real(nPUMPframe(3))/real(nframe),nCHEMframe(1),nCHEMframe(2),nCHEMframe(3),real(nCHEMframe(1))/real(nframe),real(nCHEMframe(2))/real(nframe),real(nCHEMframe(3))/real(nframe)
          nPUMPframe = 0
          nCHEMframe = 0
        endif
 
    end do

    iframe = mod(nframextc(ifile),nframe)

    nPUMPframe(2) = nPUMPframe(1)+nPUMPframe(2)
    nPUMPframe(3) = nPUMPframe(2)+nPUMPframe(3)
    nCHEMframe(2) = nCHEMframe(1)+nCHEMframe(2)
    nCHEMframe(3) = nCHEMframe(2)+nCHEMframe(3)
    write(17, '(f12.6,I,I,I,f12.6,f12.6,f12.6,I,I,I,f12.6,f12.6,f12.6)') traj % time, nPUMPframe(1),nPUMPframe(2),nPUMPframe(3), real(nPUMPframe(1))/real(nframe),real(nPUMPframe(2))/real(nframe),real(nPUMPframe(3))/real(nframe),nCHEMframe(1),nCHEMframe(2),nCHEMframe(3),real(nCHEMframe(1))/real(nframe),real(nCHEMframe(2))/real(nframe),real(nCHEMframe(3))/real(nframe)
!    write(17, '(f12.6,I,f12.6,I,f12.6)') traj % time, nPUMPframe, real(nPUMPframe)/real(iframe),nCHEMframe,real(nCHEMframe)/real(iframe)
    nPUMPframe = 0
    nCHEMframe = 0

    ! 5. Close the file
    call traj % close
    enddo

    nxtc=ixtc
    do ifile=1,nxtcfile
      write(16,'(A,A,A,I)') 'filename :',strXtcfiles(ifile),'nframe :',nframextc(ifile)
      rPUMPxtc(:) = real(nPUMPxtc(ifile,:))/real(nframextc(ifile))     
      rCHEMxtc(:) = real(nCHEMxtc(ifile,:))/real(nframextc(ifile))     
      write(16,'(A,f12.6,f12.6,A,f12.6,f12.6)') 'number of PUMP configurations :', rPUMPxtc, ' number of CHEM configurations :', rCHEMxtc
    enddo
    write(16, *) 'total number of frame         : ', nxtc
    write(16, *) 'number of PUMP configurations : ', nPUMP
    write(16, *) 'number of CHEM configurations : ', nCHEM
    write(6, *) 'total number of frame         : ', nxtc
    write(6, *) 'number of PUMP configurations : ', nPUMP
    write(6, *) 'number of CHEM configurations : ', nCHEM

    close(16) 
    close(17)

end program hbond_network
