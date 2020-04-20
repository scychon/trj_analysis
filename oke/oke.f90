!  XDR Fortran Interface xtc Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!
 
program calc_oke
 
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
  integer            :: i,j,k,l,icnt,idx,ixtc,idxmol,idxatom,idxmoltype           !counters
  integer            :: isize,ifile                  !counters
  integer            :: nummolintype                 !counters
  integer            :: ctime_0, ctime_1, ctime_rate, ctime_max          !counters_for_timer
  integer            :: narg, cptArg       !#of arg & counter of arg
  integer            :: nframe
  integer            :: nsysatoms, nmolcurr
  integer            :: idframe,ndframe,idt
  integer            :: numthread,threadid
  integer,dimension(3)   :: idxs
  integer time_array_0(8), time_array_1(8)
  real start_time, finish_time


  real*8             :: delr=0.02d0,rmax=3.0d0,boxsize=0.3d0      !
  real*8             :: dt,dt0,t1,t2,tstart,tend
  real*8             :: mass, molmass, dr, rclustcat, rclustani,charge,molcharge
  real*8             :: murotsq, mutranssq, murottrans, volume
  real*8             :: musq,cosval
  real*8             :: rotsq,transsq,rottrans,rotval,transval
  real*8             :: debye, eps0, kb, enm, multiple, qelec,mult2,multiConduct
  real*8             :: volavg, muavg
  real*8,dimension(3)    :: sysbox,dist,box,mutot,costot
  real*8,dimension(3)    :: xpos,compos   ! temporary position vector for x, center of mass(com)
  real*8,dimension(3)    :: murot, mutrans
  real*8,dimension(300)  :: rbin
  real*8, allocatable,dimension(:,:,:) :: comtraj,comtrajcat,comtrajani,mutraj,mutrajcat,mutrajani,temp   ! trajectory of com of each molecules (frame idx, mol idx, xyz)
  real*8, allocatable,dimension(:,:,:) :: costraj,costrajmol,cossqtrajmol   ! trajectory of dipole alignment vector (frame idx, mol idx, xyz)
  real*8, allocatable,dimension(:,:,:) :: comUnwrapTraj   ! unwrapped trajectory of com of each molecules (frame idx, mol idx, xyz)
  real*8, allocatable,dimension(:,:) :: comMolt1,comMolt2,comMolDiff   ! array for com of each molecules at specific time (mol idx, xyz)
  real*8, allocatable,dimension(:,:) :: boxtraj,tempbox      !matrix for # of water molecules along z direction
  real*8, allocatable,dimension(:,:) :: mutrajrot, mutrajtrans      !matrix for murot,mutran
  real*8, allocatable,dimension(:,:) :: xposframe   !matrix for positions in the frame
  real*8, allocatable,dimension(:) :: temptime     !matrix for timestamps
  real*8, allocatable,dimension(:) :: Sc     !matrix for q*xyz for each timestep
  real*8, allocatable,dimension(:) :: molTypeIdxs     ! array for molecule type idx of each moleclue ( molecule idx )
  real*8, allocatable,dimension(:) :: nmolMolType    ! array for number of molecules in each moltype ( moltype idx )

  ! matrices for time stamps
  real*8, allocatable,dimension(:) :: timestamp, timestep     !matrix for timestamps

  ! array of count of frame pairs having time difference dt
  integer,allocatable,dimension(:)   :: nDiffTime   !  (time difference dt in unit of time step)
  integer,allocatable,dimension(:,:)   :: nDiffTimeThread   !  (time difference dt in unit of time step)
  ! array of dt in unit of dt0
  integer,allocatable,dimension(:,:)   :: idtmat   !  (ixtc, jxtc)

  ! for conductivity calculation
  ! matrices
  real*8, allocatable,dimension(:,:,:) :: xqcomTraj     ! matrix for q*xyz_com of each molecule at each timestep (frame idx, mol idx, q*xyz)
  real*8, allocatable,dimension(:,:,:) :: xqAtomsTraj     ! matrix for q*xyz of each atom at each timestep (frame idx, atom idx, q*xyz)
  real*8, allocatable,dimension(:) :: xqAtomsCFTime     ! matrix for conduct correlation function of each atom at each dt (time difference dt in unit of time step, atom idx)
  real*8, allocatable,dimension(:,:) :: xqAtomsCFTimeThread     ! matrix for conduct correlation function of each atom at each dt (time difference dt in unit of time step, atom idx)
  real*8, allocatable,dimension(:) :: xqcomCFTime     ! matrix for conduct correlation function of each molecule at each dt (time difference dt in unit of time step, mol idx, q*xyz)
  real*8, allocatable,dimension(:,:) :: xqcomCFTimeThread     ! matrix for conduct correlation function of each molecule at each dt (time difference dt in unit of time step, mol idx, q*xyz)
  real*8, allocatable,dimension(:,:) :: xqcomDiff     ! matrix for difference of xq_xyz of each molecules at each time step combination (ixtc,jxtc,idxmol,3)
  real*8,dimension(3)    :: xqdiff1,xqdiff2


  ! for MSD calculation
  real*8                    :: msd, msd_axis
  real*8,dimension(3)       :: xdiff
  real*8, allocatable,dimension(:,:,:)  :: msdTime      ! translational mean square displacement of each molecule types at each time difference (time difference dt in unit of time step, x;y;z;tot, moltype idx i)
  real*8, allocatable,dimension(:,:,:,:)  :: msdTimeThread      ! translational mean square displacement of each molecule types at each time difference (time difference dt in unit of time step, x;y;z;tot, moltype idx i)

  ! for rotational ACF
  ! matrices for unit normal vector of molecules
  ! c1, c21, c22 forms the plane for imidazolium ring, b or p, f1, f2 forms the plane for bf4 or pf6
  real*8                    :: rotacf
  real*8, dimension(3,3)    :: xPlaneAtom   ! temporary position vector for the three atoms that compose the molecular plane
  real*8, allocatable,dimension(:,:,:) :: unitNormMolTraj     ! unit normal vector of each molecules at each timestep (frame idx, mol idx, vec_xyz)
  real*8, allocatable,dimension(:,:) :: unitNormMolt1, unitNormMolt2     ! unit normal vector of each molecules at each timestep (frame idx, mol idx, vec_xyz)
  !real*8, allocatable,dimension(:,:,:) :: rotacf     ! rotational ACF <ui(t2).ui(t1)> of each molecules at each time step combination (ixtc, jxtc, idxmol)
  real*8, allocatable,dimension(:,:) :: rotACFTime,rotACFTimeP2     ! rotational ACF <ui(t2).ui(t1)> of each molecule types at each time difference (time difference dt in unit of time step, moltype idx i )
  real*8, allocatable,dimension(:,:,:) :: rotACFTimeThread,rotACFTimeP2Thread     ! rotational ACF <ui(t2).ui(t1)> of each molecule types at each time difference (time difference dt in unit of time step, moltype idx i )
  real*8, allocatable,dimension(:) :: rotACFt0     ! rotational ACF <ui(t2).ui(t1)> of each molecule types at each time difference (time difference dt in unit of time step, moltype idx i )

  character(len=30)          :: strfmt, strfmtmsd
  character(len=256)         :: strout,strtitle
  character(len=256)         :: str, line
  character(len=256),allocatable,dimension(:) :: aStrArgs
! ----------------------------------------------------------------

  murotsq = 0
  mutranssq = 0
  murottrans = 0
  muavg = 0
  nframe = 0
  nsysatoms = 0
  volavg = 0
  volume = 0

  write(*,*) 'get omp num threads'
  !$OMP PARALLEL
  numthread = omp_get_num_threads()
  !$OMP END PARALLEL
  write(6,*) 'Using ',numthread,' number of threads'

  debye = 3.33564E-30
  eps0 = 8.85418781762E-12
  kb = 1.380648813E-23
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
    write(*,'(A,A)') 'topology file is read ! ',trim(strTopFile)

    nmolsys = sys % numsysmol
    nmoltype = sys % nummoltype
    allocate(molTypeIdxs(nmolsys))
    allocate(nmolMolType(nmoltype))
    idxmol = 0
    do i=1, sys % nummoltype
      nmolMolType(i) = sys % moltype(i) % nummol
      nummolintype = nmolMolType(i)
      nsysatoms = nsysatoms + nmolMolType(i) * sys % moltype(i) % numatom
      do j=1, nummolintype
        idxmol = idxmol + 1
        molTypeIdxs(idxmol) = i
      enddo
    enddo
    allocate(xposframe(3,nsysatoms))

    ! initialize all matrices for first trajectory 
    allocate(comMolt1(nmolsys,3))
    !allocate(comMolt2(nmolsys,3))
    allocate(comMolDiff(nmolsys,3))
      allocate(comtraj(nsize,nmolsys,3))
      allocate(comUnwrapTraj(nsize,nmolsys,3))
      allocate(mutraj(nsize,nmolsys,3))
      allocate(costraj(3,nmolsys,nsize))
      allocate(unitNormMolTraj(nsize,nmolsys,3))
!      allocate(xqAtomsTraj(nsize,nsysatoms,3))
      allocate(boxtraj(nsize,3))
      allocate(timestamp(nsize))
      allocate(timestep(nsize))
    write(6,*) 'comMolt1 & comMolDiff allocated' 

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
      !stop
    endif

    OPEN (UNIT=6,FORM='FORMATTED')
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
!        write(*,*) 'volume', volume, trajin % box(1,1), trajin % box(2,2), trajin % box(3,3)

        ! check the size of matrices 
        isize = size(comtraj(:,1,1))
        if(isize .lt. ixtc) then
          write(6,*) 'trajectory is longer than ', nsize
          call expand3D(comtraj,nsize,0,0)
          call expand3D(comUnwrapTraj,nsize,0,0)
          call expand3D(mutraj,nsize,0,0)
          call expand3D(costraj,0,0,nsize)
          call expand3D(unitNormMolTraj,nsize,0,0)
        !  call expand3D(xqAtomsTraj,nsize,0,0)
          call expand2D(boxtraj,nsize,0)
          call expand1D(timestamp,nsize)
          call expand1D(timestep,nsize)
          nsize = nsize*2
        endif

        ! record current time
        timestamp(ixtc) = trajin % time
        timestep(ixtc) = trajin % STEP


        do k=1,3
          sysbox(k) = trajin % box(k,k)
          boxtraj(ixtc,k) = sysbox(k)
        enddo

        ! Compute and store the trajectory of the center of mass postion,
        ! the rotational principle vector, and the dipole moment vector
        idxatom = 0
        idxmol = 0
        xposframe = trajin % pos
        do i=1, sys % nummoltype
          !$OMP PARALLEL &
          !$OMP   DEFAULT (FIRSTPRIVATE) &
          !$OMP   SHARED (unitNormMolTraj, comtraj, mutraj,costraj)
          !$OMP DO
          do j=1, sys % moltype(i) % nummol
            compos = 0
            mutot = 0
            idxs = sys % moltype(i) % molPlaneAtomIdxs

            do k=1, sys % moltype(i) % numatom
              idx = idxatom + (j-1) * sys % moltype(i) % numatom + k
              mass = sys % moltype(i) % atommass(k)
              charge = sys % moltype(i) % atomcharge(k)
              xpos(:) = xposframe(:,idx)
              compos(:) = compos(:) + xpos(:)*mass
            !  xqAtomsTraj(ixtc,idxatom,:) = xpos(:)*charge
              mutot(:) = mutot(:) + xpos(:)*charge
              if (k .eq. idxs(1) ) then
                  xPlaneAtom(1,:) = xpos(:)
              elseif (k .eq. idxs(2) ) then
                  xPlaneAtom(2,:) = xpos(:)
              elseif (k .eq. idxs(3) ) then
                  xPlaneAtom(3,:) = xpos(:)
              end if
            enddo
            unitNormMolTraj(ixtc,idxmol+j,:) = vecUnitNorm( xPlaneAtom(1,:)-xPlaneAtom(2,:), xPlaneAtom(1,:)-xPlaneAtom(3,:) )

            molmass = sys % moltype(i) % molarmass
            compos(:) = compos(:)/molmass
            comtraj(ixtc,idxmol+j,:) = compos(:)
            mutraj(ixtc,idxmol+j,:) = mutot(:)
            costraj(:,idxmol+j,ixtc) = mutot(:)/norm(mutot)
!            write(*,*) comtraj(ixtc,idxmol+j,:)
          enddo
          !$OMP END DO
          !$OMP END PARALLEL
          idxmol = idxmol + sys % moltype(i) % nummol
          idxatom = idxatom + sys % moltype(i) % nummol * sys % moltype(i) % numatom
        enddo 

        ! unwrap the trajectory
        if(ixtc .eq. 1) then
          comUnwrapTraj(ixtc,:,:) = comtraj(ixtc,:,:)
        else if (ixtc .gt. 1) then
          !$OMP PARALLEL &
          !$OMP   DEFAULT (FIRSTPRIVATE) &
          !$OMP   SHARED (comtraj, comMolDiff, comUnwrapTraj)
          !$OMP DO
          do i=1, nmolsys
            comMolDiff(i,:) = comtraj(ixtc,i,:) - comtraj(ixtc-1,i,:)
            comMolDiff(i,1) = comMolDiff(i,1) - NINT(comMolDiff(i,1)/sysbox(1))*sysbox(1)
            comMolDiff(i,2) = comMolDiff(i,2) - NINT(comMolDiff(i,2)/sysbox(2))*sysbox(2)
            comMolDiff(i,3) = comMolDiff(i,3) - NINT(comMolDiff(i,3)/sysbox(3))*sysbox(3)
            comUnwrapTraj(ixtc,i,:) = comUnwrapTraj(ixtc-1,i,:)+comMolDiff(i,:)
          enddo
          !$OMP END DO
          !$OMP END PARALLEL
        endif 

        ! check number of atoms in system
        if(nsysatoms .ne. idxatom) then
          write(*,*) 'number of atoms in ',ixtc,'th trajectory does not match other frames'
          write(*,*) 'nSysAtoms = ',nsysatoms,'  idxatom = ', idxatom
          nsysatoms = idxatom
          !stop
        endif

        write(6,100,advance='no') achar(13), ixtc,'th frame has finished  '
 100    FORMAT(A,I8,A)

        call trajin % read
 
    end do

    ! 5. Close the file
    call trajin % close
    end do

    nxtc = ixtc
    nframe = nframe + nxtc
    volavg = volavg * 1e-27 / nxtc
    multiple = multiple / (temperature * volavg * 3)
    multiConduct = qelec**2 / (6. * kb * temperature * volavg)
    write(*,*) 'temp mult vol', multiConduct, temperature, volavg

    write(6,*) 'start analyzing the data' 

    ! check the size of trajectory matrix
    isize = size(comtraj(:,1,1))
    if(isize .gt. nxtc) then
      call shrink3D(comtraj,nxtc,nmolsys,3)
      call shrink3D(comUnwrapTraj,nxtc,nmolsys,3)
      call shrink3D(mutraj,nxtc,nmolsys,3)
      call shrink3D(costraj,3,nmolsys,nxtc)
      call shrink3D(unitNormMolTraj,nxtc,nmolsys,3)
    !  call shrink3D(xqAtomsTraj,nxtc,nsysatoms,3)
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

    write(6,*) 'dt is ', dt
    write(6,*) 'dt0 is ', dt0

    write(6,*) 'ndframe is set' , ndframe
    ! generate and initialize time difference matrices
!    allocate(idtmat(nxtc-1,nxtc))
    write(6,*) 'idtmat allocated' 
    allocate(unitNormMolt1(nmolsys,3))
    allocate(unitNormMolt2(nmolsys,3))
    write(6,*) 'unitNormMolt1 & unitNormMolDiff allocated' 

!    allocate(xqAtomsCFTime(ndframe))
    allocate(xqcomCFTime(ndframe))
    write(6,*) 'xqcomCFTime allocated' 
    allocate(xqcomDiff(nmolsys,3))
    write(6,*) 'xqcomDiff allocated' 
    allocate(msdTime(ndframe,5,nmoltype+1))
    write(6,*) 'msdTime allocated' 
    allocate(rotACFTime(ndframe,nmoltype+1))
    write(6,*) 'rotACFTime allocated' 
    allocate(rotACFTimeP2(ndframe,nmoltype+1))
    write(6,*) 'rotACFTimeP2 allocated' 
 !   allocate(rotacf(nxtc-1,nxtc,nmolsys))
!    allocate(xqAtomsDiff(nxtc-1,nxtc,nsysatoms,3))
    allocate(xqcomTraj(nxtc,nmolsys,3))
    write(6,*) 'xqcomTraj allocated' 
 !   allocate(nACFTime(ndframe))
    allocate(nDiffTime(ndframe))
    write(6,*) 'nDiffTime allocated' 
    allocate(rotACFt0(nmoltype+1))
    allocate(costrajmol(3,nmoltype+1,nxtc))
    allocate(cossqtrajmol(3,nmoltype+1,nxtc))
    write(6,*) 'costraj allocated' 
!    xqAtomsCFTime = 0
    xqcomCFTime = 0
    xqcomDiff = 0
!    xqAtomsDiff = 0
    msdTime = 0
    rotACFTime = 0
    rotACFTimeP2 = 0
    rotACFt0 = 1
 !   rotacf = 0
    xqcomTraj = 0
 !   nACFTime = 0
    nDiffTime = 0

    costrajmol=0
    cossqtrajmol=0

    write(6,*) 'generating averaged molecular dipole data for each molecule types' 
    idxmol = 0
    !$OMP PARALLEL &
    !$OMP   DEFAULT (FIRSTPRIVATE) &
    !$OMP   SHARED (costraj,costrajmol,cossqtrajmol)
    !$OMP DO
    do ixtc=1,nxtc
      idxmol = 0
      do i=1, sys % nummoltype
        do j=1, sys % moltype(i) % nummol
          costot(:) = costraj(:,idxmol+j,ixtc)
          costrajmol(:,i,ixtc) = costrajmol(:,i,ixtc) + costot(:)
          cossqtrajmol(:,i,ixtc) = cossqtrajmol(:,i,ixtc) + costot(:)*costot(:)
        enddo
        nummolintype = sys % moltype(i) % nummol  
        idxmol = idxmol + nummolintype
        costrajmol(:,nmoltype+1,ixtc) = costrajmol(:,nmoltype+1,ixtc) + costrajmol(:,i,ixtc)
        cossqtrajmol(:,nmoltype+1,ixtc)=cossqtrajmol(:,nmoltype+1,ixtc)+cossqtrajmol(:,i,ixtc)
        costrajmol(:,i,ixtc) = costrajmol(:,i,ixtc)/nummolintype
        cossqtrajmol(:,i,ixtc) = cossqtrajmol(:,i,ixtc)/nummolintype
      enddo
      costrajmol(:,nmoltype+1,ixtc) = costrajmol(:,nmoltype+1,ixtc)/idxmol
      cossqtrajmol(:,nmoltype+1,ixtc) = cossqtrajmol(:,nmoltype+1,ixtc)/idxmol
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call backupfile(strCOSFile)
    call backupfile(strCOSSqFile)
    open(18,file=strCOSFile)
    open(19,file=strCOSSqFile)
    write(strfmt,'("( F12.3, ",I0,"ES15.7 )")') (nmoltype+1)*3
    write(18,'(A)') '#time,    cos_x,     cos_y,     cos_z'
    write(19,'(A)') '#time,    cos_x,     cos_y,     cos_z'

    do ixtc=1,nxtc
        write(18,strfmt) dt0*ixtc,costrajmol(:,:,ixtc)
        write(19,strfmt) dt0*ixtc,cossqtrajmol(:,:,ixtc)
    enddo

    write(6,*) 'finished writing output files'
    close(6)
    close(18)
    close(19)

end program calc_oke
