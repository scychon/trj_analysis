!  XDR Fortran Interface xtc Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!
 
program calc_xqCF_rotACF
 
    ! 1. Use the xdr interface
    use traj, only: trajfile
    use topol, only: topfile
    use fvector
    use filehandle
    use variables
    use omp_lib
 
    implicit none

    ! 2. Declare a variable of type xtcfile
    type(trajfile) :: trajin
    type(topfile) :: sys

!  --------------------------------------------------------------------------
  integer            :: nbins=200                   !# of bins 
  integer            :: nskip=100                   !# of bins 
  integer            :: nsize=10000                 !# of bins 
  integer            :: i,j,k,l,icnt,idx,ixtc,idxmol,idxatom,idxatomtype,idxmoltype,idxdrude           !counters
  integer            :: isize,ifile                  !counters
  integer            :: narg, cptArg       !#of arg & counter of arg
  integer            :: nframe
  integer            :: numthread
  integer            :: nsysatoms, nmolcurr, natomcurr, ndrudepaircurr, natomtypes
  integer            :: drudepairidx
  integer            :: ndframe,idt
  integer            :: idxstart,idxend,idxskip
  integer time_array_0(8), time_array_1(8), time_array_t(8)


  real*8             :: infinity = 1e30
  real*8             :: dt,dt0
  real*8             :: tempval
  real*8             :: mass, molmass, charge, drudemass
  real*8             :: murotsq, mutranssq, murottrans, volume
  real*8             :: debye, eps0, kb, R, mult, enm, multiple, qelec,mult2,multiConduct
  real*8             :: volavg, muavg
  real*8,dimension(3)    :: xpos,xdrudepos,compos,drudepairpos   ! temporary position vector for x, center of mass(com)
  real*8, allocatable,dimension(:,:) :: comtraj, comtrajmol   ! trajectory of com of each molecules (frame idx, mol idx, xyz)
  real*8, allocatable,dimension(:,:) :: atomveltraj   ! trajectory of com of each molecules (frame idx, mol idx, xyz)
  real*8, allocatable,dimension(:,:) :: boxtraj      !matrix for # of water molecules along z direction
  real*8, allocatable,dimension(:,:) :: xposframe   !matrix for positions in the frame
  real*8, allocatable,dimension(:) :: keframe     !matrix for kinetic energy of each atom in the frame
  real*8, allocatable,dimension(:) :: molTypeIdxs     ! array for molecule type idx of each moleclue ( molecule idx )
  real*8, allocatable,dimension(:) :: nmolMolType    ! array for number of molecules in each moltype ( moltype idx )

  ! matrices for time stamps
  real*8, allocatable,dimension(:) :: timestamp, timestep     !matrix for timestamps



  character(len=30)          :: strfmt
  character(len=256)         :: str
  character(len=256),allocatable,dimension(:) :: aStrArgs
! ----------------------------------------------------------------

  l = 0
  icnt = 0
  idxmoltype = 0
  nbins = 200
  nskip = 100
  murotsq = 0
  mutranssq = 0
  murottrans = 0
  muavg = 0
  nframe = 0
  nsysatoms = 0
  natomtypes = 0
  volavg = 0
  volume = 0

  write(*,*) 'get omp num threads'
  numthread = 0
  !$OMP PARALLEL
  numthread = omp_get_num_threads()
  !$OMP END PARALLEL
  write(6,*) 'Using ',numthread,' number of threads'

  debye = 3.33564E-30
  eps0 = 8.85418781762E-12
  kb = 1.380648813E-23
  R = 8.3144598
  enm = 0.020819434
  qelec = 1.60217656535E-19

  multiple = debye**2/(enm**2 * eps0 * kb)
  mult2 = qelec**2 * 1e-18/(eps0*kb)
  multiConduct = qelec**2 / kb

  write(*,*) multiple, mult2,multiConduct

  !Check if any arguments are found
  narg=command_argument_count()
  !Loop over the arguments
  if(narg>0)then
    !loop across options
    allocate(aStrArgs(narg))
    do cptArg=1,narg
      call get_command_argument(cptArg,str)
      aStrArgs(cptArg) = str
    enddo
    !assign in and out files
    call getargs(aStrArgs)
  else
    write(*,*) 'no arguments are given'
    write(*,*) 'usage : calc_cond -f param.dat -o outfile.dat'
    stop
  endif


    ! Read param file to locate strXtcfiles
    call readparam

    ! 3. Read the topology file and initialize atomic mass data for system
    write(*,'(A,A)') 'reading topoly file ',trim(strTopFile)
    write(*,*) ''
    call sys % init(strTopFile)

    nmolsys = sys % numsysmol
    nmoltype = sys % nummoltype
    allocate(molTypeIdxs(nmolsys))
    allocate(nmolMolType(nmoltype))
    idxmol = 0
    do i=1, sys % nummoltype
      nmolMolType(i) = sys % moltype(i) % nummol
      nsysatoms = nsysatoms + nmolMolType(i) * sys % moltype(i) % numatom
      natomtypes = natomtypes + sys % moltype(i) % numatom
      do j=1, nmolMolType(i)
        idxmol = idxmol + 1
        molTypeIdxs(idxmol) = i
      enddo
      write(6,*) 'moltype i, nmol', i, nmolMolType(i)
    enddo
    allocate(xposframe(3,nsysatoms))
    allocate(keframe(nsysatoms))

    ! initialize all matrices for first trajectory 
      allocate(comtraj(nsize,nmolsys))
      allocate(comtrajmol(nsize,nmoltype))
      allocate(atomveltraj(nsize,natomtypes))
      allocate(boxtraj(nsize,3))
      allocate(timestamp(nsize))
      allocate(timestep(nsize))

    ! backup temperature file and create new one
    call backupfile(strTempFile)
    open(18,file=strTempFile)

    call date_and_time(values=time_array_0)
    write(6,*) 'Start time : ', time_array_0
    ! 3. Initialize it with the name of trajin file you want to read in.
    ixtc=0
    do ifile=1,nxtcfile
    strInFile = trim(strXtcfiles(ifile))
    write(*,'(A,A)') 'reading trajectory file ',trim(strInFile)
    write(*,*) ''
    call trajin % init(strInFile)
 
    ! 4. Read in each configuration. Everything is stored in the xtcfile
    ! type (precision, time,
    !    step, no of atoms, positions, etc.). Look in the xtc module for
    !    more details.
    !    You can save the positions in the loop for your calculations in
    !    another array, or 
    !    do your calculations after each read.
 
    call trajin % read

    ! check number of atoms in system
    if(nsysatoms .ne. trajin % NATOMS) then
      write(*,*) 'number of atoms in ',ifile,'th trajectory file does not match other trajectories'
      write(*,*) 'nSysAtoms = ',nsysatoms,'  trajin % NATOMS = ', trajin % NATOMS
      stop
    endif

    write(*,*) 'here is the beginning point ! '
    OPEN (UNIT=6,FORM='FORMATTED',CARRIAGECONTROL='FORTRAN')
        ! Just an example to show what was read in
        write(6,'(a,f12.6,a,i0)') " Time (ps): ", trajin % time, "  Step: ", trajin % STEP
        write(6,'(a,f12.6,a,i0)') " Precision: ", trajin % prec, "  No. Atoms: ", trajin % NATOMS
        ! This is the same order as found in the GRO format fyi
        write(6,'(9f9.5)')  trajin % box(1,1), trajin % box(2,2), trajin % box(3,3), &
                            trajin % box(1,2), trajin % box(1,3), & 
                            trajin % box(2,1), trajin % box(2,3), &
                            trajin % box(3,1), trajin % box(3,2)
        write(6,*) '' 
    do while ( trajin % STAT == 0 )
        ixtc = ixtc + 1
        volume = trajin % box(1,1) * trajin % box(2,2) * trajin % box(3,3)
        volavg = volavg + volume
        !write(*,*) 'volume', volume, trajin % box(1,1), trajin % box(2,2), trajin % box(3,3)

        ! check the size of matrices 
        isize = size(comtraj(:,1))
        if(isize .lt. ixtc) then
          write(6,*) 'trajectory is longer than ', nsize
          call expand2D(comtraj,nsize,0)
          call expand2D(comtrajmol,nsize,0)
          call expand2D(atomveltraj,nsize,0)
          call expand2D(boxtraj,nsize,0)
          call expand1D(timestamp,nsize)
          call expand1D(timestep,nsize)
          nsize = nsize*2
        endif

        ! record current time
        timestamp(ixtc) = trajin % time
        timestep(ixtc) = trajin % STEP

        !write(*,*) timestamp(ixtc)
        !write(*,*) 'initial time recoded !'

        do k=1,3
          boxtraj(ixtc,k) = trajin % box(k,k)
        enddo

        idxatom = 0
        idxmol = 0
        idxatomtype = 0
        xposframe = trajin % pos
        !write(6,*) 'position set !'
        do i=1, sys % nummoltype
          nmolcurr = sys % moltype(i) % nummol
          natomcurr = sys % moltype(i) % numatom
          ndrudepaircurr = sys % moltype(i) % numdrudepair
          molmass = sys % moltype(i) % molarmass
          !write(6,*) 'nummol', sys % moltype(i) % nummol
          !$OMP PARALLEL &
          !$OMP   DEFAULT (FIRSTPRIVATE) &
          !$OMP   SHARED (comtraj,keframe)
          !$OMP DO
          do j=1, nmolcurr
            compos = 0
            idxdrude = 0
            do k=1, natomcurr
              idx = idxatom + (j-1) * natomcurr + k
              mass = sys % moltype(i) % atommass(k)
              charge = sys % moltype(i) % atomcharge(k)
              xpos(:) = xposframe(:,idx)
              !write(6,*) xpos
              compos(:) = compos(:) + xpos(:)*mass
              drudepairidx = sys % moltype(i) % atomdrudepairidx(k)
              if (drudepairidx .gt. k) then
                xdrudepos(:) = xposframe(:,idx-k+drudepairidx)
                drudemass = sys % moltype(i) % atommass(drudepairidx)
                drudepairpos(:) = mass * xpos(:) + drudemass * xdrudepos(:)
                keframe(idx) = dot_product(drudepairpos,drudepairpos)/(mass+drudemass)
              else if (drudepairidx .gt. 0) then
                xdrudepos(:) = xposframe(:,idx-k+drudepairidx)
                drudemass = sys % moltype(i) % atommass(drudepairidx)
                drudepairpos(:) = mass * xpos(:) + drudemass * xdrudepos(:)
                keframe(idx) = mass * dot_product(xpos,xpos) + drudemass * dot_product(xdrudepos,xdrudepos)
                keframe(idx) = keframe(idx) - dot_product(drudepairpos,drudepairpos)/(mass+drudemass)
              else
                keframe(idx) = mass * dot_product(xpos,xpos)
              endif
            enddo
            comtraj(ixtc,idxmol+j) = dot_product(compos,compos)/molmass
          !  write(6,*) j,comtraj(ixtc, idxmol+j)
          enddo
          !$OMP END DO
          !$OMP END PARALLEL

          comtrajmol(ixtc,i) = sum(comtraj(ixtc,idxmol+1:(idxmol+nmolcurr)))/nmolcurr
          if (comtrajmol(ixtc,i) .gt. infinity) then
            write(6,*) 'comtrajmol is infinity !'
            write(6,*) comtraj(ixtc,idxmol+1:(idxmol+nmolcurr))
            write(6,*) ''
          endif

          idxend = idxatom + nmolcurr*natomcurr
          idxskip = natomcurr
          !write(*,*) 'comtraj read !'
          !!$OMP PARALLEL &
          !!$OMP   DEFAULT (FIRSTPRIVATE) &
          !!$OMP   SHARED (atomveltraj)
          !!$OMP DO
          do j=1, sys % moltype(i) % numatom
            !write(*,*) 'aomtype', j, keframe((idxatom+j):(idxatom+nmolcurr*natomcurr):natomcurr)
            idxstart = idxatom + j
            tempval = sum(keframe(idxstart:idxend:idxskip))/float(nmolcurr)
            !tempval = sum(keframe((idxatom+j):(idxatom+nmolcurr*natomcurr):natomcurr))/nmolcurr
            atomveltraj(ixtc,idxatomtype+j) = tempval
           ! write(*,*) atomveltraj(ixtc,idxatomtype+j)
          enddo
          !!$OMP END DO
          !!$OMP END PARALLEL
          !write(*,*) 'atomveltraj read !'
          idxmol = idxmol + sys % moltype(i) % nummol
          idxatom = idxatom + sys % moltype(i) % nummol * sys % moltype(i) % numatom
          idxatomtype = idxatomtype + sys % moltype(i) % numatom
        enddo 

        !write(*,*) 'frame data read !'

        ! check number of atoms in system
        if(nsysatoms .ne. idxatom) then
          write(*,*) 'number of atoms in ',ixtc,'th trajectory does not match other frames'
          write(*,*) 'nSysAtoms = ',nsysatoms,'  idxatom = ', idxatom
          stop
        endif

        write(6,100) ixtc,'th frame has finished  ' 
 100    FORMAT('+', I5,A)  

        call trajin % read
 
    end do

    ! 5. Close the file
    call trajin % close
    end do
    call date_and_time(values=time_array_1)
    write(6,*) 'End time : ', time_array_1

    nxtc = ixtc
    nframe = nframe + nxtc
    volavg = volavg * 1e-27 / nxtc
    multiple = multiple / (temperature * volavg * 3)
    mult = real(1e3) / (3.d0 * R)
    multiConduct = qelec**2 / (6. * kb * temperature * volavg)
    write(*,*) 'temp mult vol', temperature, mult, volavg

    write(6,*) 'start analyzing the data' 

    ! check the size of trajectory matrix
    isize = size(comtraj(:,1))
    if(isize .gt. nxtc) then
      call shrink2D(comtraj,nxtc,nmolsys)
      call shrink2D(comtrajmol,nxtc,nmoltype)
      call shrink2D(atomveltraj,nxtc,natomtypes)
      call shrink2D(boxtraj,nxtc,3)
      call shrink1D(timestamp,nxtc)
      call shrink1D(timestep,nxtc)
    endif

    write(6,*) 'matrix size has adjusted' 

    ! calculate size of time difference matrices
    ! assume the first two steps have minimum dt
    ! and the trajectories are time-ordered
    dt0 = timestamp(2)-timestamp(1)
    dt = timestamp(nxtc)-timestamp(1)
    ndframe = int(dt/dt0 + 0.00001)

    write(6,*) 'ndframe is set' , ndframe

    write(6,*) 'data analysis done'
    

    write(6,*) 'start generating output files'


    write(strfmt,'("( F12.3, ",I0,"ES15.7 )")') nmoltype + natomtypes
    do idt=1,ndframe
      write(18,strfmt) dt0*idt,comtrajmol(idt,:)*mult,atomveltraj(idt,:)*mult
    enddo


    write(6,*) 'finished writing output files'
    close(6)
    close(18)

end program calc_xqCF_rotACF
