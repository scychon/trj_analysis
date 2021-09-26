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
  integer            :: i,j,k,l,icnt,idx,ixtc,idxmol,idxatom,idxmoltype,idxcurr,idxcurr2      !counters
  integer            :: idxlast                     ! last index for xqcomCFTime
  integer            :: isize,ifile                  !counters
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
  real*8             :: rotsq,transsq,rottrans,rotval,transval
  real*8             :: debye, eps0, kb, enm, multiple, qelec,mult2,multiConduct
  real*8             :: volavg, muavg
  real*8,dimension(3)    :: sysbox,dist,box,mutot
  real*8,dimension(3)    :: xpos,compos   ! temporary position vector for x, center of mass(com)
  real*8,dimension(3)    :: murot, mutrans
  real*8,dimension(300)  :: rbin
  real*8, allocatable,dimension(:,:,:) :: comtraj,comtrajcat,comtrajani,mutraj,mutrajcat,mutrajani,temp   ! trajectory of com of each molecules (frame idx, mol idx, xyz)
  real*8, allocatable,dimension(:,:,:) :: comUnwrapTraj   ! unwrapped trajectory of com of each molecules (frame idx, mol idx, xyz)
  real*8, allocatable,dimension(:,:) :: comMolt1,comMolt2,comMolDiff   ! array for com of each molecules at specific time (mol idx, xyz)
  real*8, allocatable,dimension(:,:) :: boxtraj,tempbox      !matrix for # of water molecules along z direction
  real*8, allocatable,dimension(:,:) :: mutrajrot, mutrajtrans      !matrix for murot,mutran
  real*8, allocatable,dimension(:,:) :: xposframe   !matrix for positions in the frame
  real*8, allocatable,dimension(:) :: temptime     !matrix for timestamps
  real*8, allocatable,dimension(:) :: Sc     !matrix for q*xyz for each timestep
  integer, allocatable,dimension(:) :: molTypeIdxs     ! array for molecule type idx of each moleclue ( molecule idx )
  integer, allocatable,dimension(:) :: nmolMolType    ! array for number of molecules in each moltype ( moltype idx )

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
  real*8, allocatable,dimension(:,:) :: xqcomCFTime     ! matrix for conduct correlation function of each molecule at each dt (time difference dt in unit of time step, mol idx, q*xyz)
  real*8, allocatable,dimension(:,:,:) :: xqcomCFTimeThread     ! matrix for conduct correlation function of each molecule at each dt (time difference dt in unit of time step, mol idx, q*xyz)
  real*8, allocatable,dimension(:,:) :: xqcomDiff     ! matrix for difference of xq_xyz of each molecules at each time step combination (ixtc,jxtc,idxmol,3)
  real*8, allocatable,dimension(:) :: qmolMolType     ! matrix for molecular charge of each moltype (nmoltype)
  real*8,dimension(3)    :: xqdiff1,xqdiff2,xqdiffsum


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

  character(len=30)          :: strfmt, strfmtmsd, strfmtconduct
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
    allocate(qmolMolType(nmoltype))
    idxmol = 0
    do i=1, sys % nummoltype
      nmolMolType(i) = sys % moltype(i) % nummol
      nsysatoms = nsysatoms + nmolMolType(i) * sys % moltype(i) % numatom
      nmolcurr = nmolMolType(i)
      do j=1, nmolcurr
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
!          !$OMP PARALLEL &
!          !$OMP   DEFAULT (FIRSTPRIVATE) &
!          !$OMP   SHARED (unitNormMolTraj, comtraj, mutraj)
!          !$OMP DO
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
!            write(*,*) comtraj(ixtc,idxmol+j,:)
          enddo
!          !$OMP END DO
!          !$OMP END PARALLEL
          idxmol = idxmol + sys % moltype(i) % nummol
          idxatom = idxatom + sys % moltype(i) % nummol * sys % moltype(i) % numatom
        enddo 

        ! unwrap the trajectory
        if(ixtc .eq. 1) then
          comUnwrapTraj(ixtc,:,:) = comtraj(ixtc,:,:)
        else if (ixtc .gt. 1) then
!          !$OMP PARALLEL &
!          !$OMP   DEFAULT (FIRSTPRIVATE) &
!          !$OMP   SHARED (comtraj, comMolDiff, comUnwrapTraj)
!          !$OMP DO
          do i=1, nmolsys
            comMolDiff(i,:) = comtraj(ixtc,i,:) - comtraj(ixtc-1,i,:)
            comMolDiff(i,1) = comMolDiff(i,1) - NINT(comMolDiff(i,1)/sysbox(1))*sysbox(1)
            comMolDiff(i,2) = comMolDiff(i,2) - NINT(comMolDiff(i,2)/sysbox(2))*sysbox(2)
            comMolDiff(i,3) = comMolDiff(i,3) - NINT(comMolDiff(i,3)/sysbox(3))*sysbox(3)
            comUnwrapTraj(ixtc,i,:) = comUnwrapTraj(ixtc-1,i,:)+comMolDiff(i,:)
          enddo
!          !$OMP END DO
!          !$OMP END PARALLEL
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
    if (bCalcTransfer) then
        idxlast = (nmoltype)*(nmoltype)+1
    else
        idxlast = 1
    endif
    allocate(xqcomCFTime(ndframe,idxlast))
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

    write(6,*) 'generating comUnwrapTraj and xqcomTraj matrix' 
    idxmol = 0
    do i=1, sys % nummoltype
      molcharge = sys % moltype(i) % molcharge
      qmolMolType(i) = molcharge
      !$OMP PARALLEL &
      !$OMP   DEFAULT (FIRSTPRIVATE) &
      !$OMP   SHARED (comtraj,comUnwrapTraj,xqcomTraj)
      !$OMP DO
      do j=1, sys % moltype(i) % nummol
        !xqcomTraj(:,idxmol+j,:) = comtraj(:,idxmol+j,:)*molcharge
        xqcomTraj(:,idxmol+j,:) = comUnwrapTraj(:,idxmol+j,:)*molcharge
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      idxmol = idxmol + sys % moltype(i) % nummol
    enddo


    write(6,*) 'start generating rotacf and xqdiff matrices'
    write(6,*) 'first frame start'
    
    call date_and_time(values=time_array_0)
    write(6,*) 'Start time : ', time_array_0
    do i=1,nxtc-1
      if ((i .gt. 100) .and. (mod(i,100) .ne. 0)) then
        cycle
      endif
      t1 = timestamp(i)
      unitNormMolt1(:,:) = unitNormMolTraj(i,:,:)
!      write(6,*) 'unitNormMol set'
!      write(6,*) 'test line'
      !$OMP PARALLEL &
      !$OMP PRIVATE(t2,dt,idt,comMolDiff,unitNormMolt2,idxmol,xqcomDiff, &
      !$OMP         xqdiff1,k,idxmoltype,nmolcurr,xdiff,msd,msd_axis,rotacf, &
      !$OMP         xqdiff2, xqdiffsum, l,idxcurr)
      !$OMP DO
      do j=i+1,nxtc
        t2 = timestamp(j)
        dt = t2 - t1
        idt = int((dt+0.00001)/dt0 + 0.0001)
        if(idt .le. 0) then
          write(6,*) 'idt is less than 0', idt, t1, t2, dt, dt0
          write(6,*) 'error might occur'
          write(6,*) 'dt/dt0', dt/dt0
          write(6,*) 'error might occur'
          write(6,*) 'error might occur'
          continue
        elseif((idt .gt. 1000) .and. (mod(idt,10) .ne. 0) ) then
          cycle
        elseif((idt .gt. 10000) .and. (mod(idt,100) .ne. 0) ) then
          cycle
        endif
        !comMolDiff(:,:) = comtraj(i,:,:)- comtraj(j,:,:)
        comMolDiff(:,:) = comUnwrapTraj(i,:,:)- comUnwrapTraj(j,:,:)
        unitNormMolt2(:,:) = unitNormMolTraj(j,:,:)

     !   nACFTime(idt) = nACFTime(idt) + 1
        nDiffTime(idt) = nDiffTime(idt) + 1

        idxmol = 0
        xqcomDiff(:,:) = xqcomTraj(i,:,:)-xqcomTraj(j,:,:)
        do k=1,nmolsys
          idxmoltype = molTypeIdxs(k)
          nmolcurr = nmolMolType(idxmoltype)
          xdiff(:) = comMolDiff(k,:)
          msd = dot_product(xdiff, xdiff)
          msd_axis = dot_product(xdiff, msd_vec)
          msd_axis = msd_axis**2
          xdiff = vecSquare(xdiff)
          msdTime(idt,1:3,idxmoltype) = msdTime(idt,1:3,idxmoltype) + xdiff(:)/nmolcurr
          msdTime(idt,1:3,nmoltype+1) = msdTime(idt,1:3,nmoltype+1) + xdiff(:)/nmolsys
          msdTime(idt,4,idxmoltype)   = msdTime(idt,4,idxmoltype) + msd_axis/nmolcurr
          msdTime(idt,4,nmoltype+1)   = msdTime(idt,4,nmoltype+1) + msd_axis/nmolsys
          msdTime(idt,5,idxmoltype)   = msdTime(idt,5,idxmoltype) + msd/nmolcurr
          msdTime(idt,5,nmoltype+1)   = msdTime(idt,5,nmoltype+1) + msd/nmolsys
          
          rotacf = dot_product(unitNormMolt1(k,:),unitNormMolt2(k,:))
          rotACFTime(idt,idxmoltype) = rotACFTime(idt,idxmoltype) + rotacf/nmolcurr
          rotACFTime(idt,nmoltype+1) = rotACFTime(idt,nmoltype+1) + rotacf/nmolsys
          rotacf = rotacf*rotacf
          rotACFTimeP2(idt,idxmoltype) = rotACFTimeP2(idt,idxmoltype) + rotacf/nmolcurr
          rotACFTimeP2(idt,nmoltype+1) = rotACFTimeP2(idt,nmoltype+1) + rotacf/nmolsys
        enddo

        if (bCalcTransfer) then
            do k=1,nmolsys
              if (qmolMoltype(molTypeIdxs(k)) .eq. 0) cycle
              xqdiff1(:) = xqcomDiff(k,:)
              do l=k,nmolsys
                if (qmolMoltype(molTypeIdxs(l)) .eq. 0) cycle
                xqdiff2(:) = xqcomDiff(l,:)
                idxcurr = molTypeIdxs(l) + (molTypeIdxs(k)-1)*(nmoltype)
                if (l==k) then
                    xqcomCFTime(idt,idxcurr) = xqcomCFTime(idt,idxcurr) + dot_product(xqdiff1,xqdiff2)
                else
                    xqcomCFTime(idt,idxcurr) = xqcomCFTime(idt,idxcurr) + 2*dot_product(xqdiff1,xqdiff2)
                endif
              enddo
            enddo
        endif
        xqdiff1 = 0
        do k=1,nmolsys
            xqdiff1(:) = xqdiff1(:) + xqcomDiff(k,:)
        enddo
        xqcomCFTime(idt,idxlast) = xqcomCFTime(idt,idxlast) + dot_product(xqdiff1,xqdiff1)

      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      write(6,100,advance='no') achar(13), i,'th frame has finished  ' 
    enddo

    call date_and_time(values=time_array_1)
    write(6,*) 'End time : ', time_array_1

      start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 &
           + time_array_0 (7) + 0.001 * time_array_0 (8)

      finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60 &
           + time_array_1 (7) + 0.001 * time_array_1 (8)
    write(6,*) 'Elapsed time : ', finish_time - start_time, ' s'

    write(6,*) 'rotacf and xqdiff matrices generated'
    write(6,*) 'generating correlation function matrices'


!    open(21,file='testout')
!    do i=1,nmolsys
!      xqdiff1(:) = xqcomTraj(1,i,:)-xqcomTraj(11,i,:)
!      write(21,*) i, xqcomTraj(1,i,:),xqcomTraj(11,i,:),dot_product(xqdiff1,xqdiff1)
!      write(21,*) i,unitNormMolTraj(1,i,:),unitNormMolTraj(2,i,:), dot_product(unitNormMolTraj(1,i,:),unitNormMolTraj(2,i,:))
!    enddo

    write(6,*) 'data analysis done'
    

    write(6,*) 'start generating output files'

    call backupfile(strMSDFile)
    call backupfile(strConductFile)
    call backupfile(strRotACFFile)
    call backupfile(strRotACFP2File)
    open(18,file=strMSDFile)
    open(19,file=strConductFile)
    open(20,file=strRotACFFile)
    open(21,file=strRotACFP2File)

!    xqAtomsCFTime = xqAtomsCFTime*multiConduct
    xqcomCFTime = xqcomCFTime*multiConduct
    if (bCalcTransfer) then
      do k=1,nmoltype-1
        do l=k+1,nmoltype
            idxcurr = l + (k-1)*nmoltype
            idxcurr2 = k + (l-1)*nmoltype
            xqcomCFTime(:,idxcurr2) = xqcomCFTime(:,idxcurr)
        enddo
      enddo
    endif


    write(18,'(A5,3E15.7)') 'tot  ',rotsq,transsq,rottrans
    write(strfmt,'("( F12.3, ",I0,"ES15.7 )")') nmoltype+1
    write(strfmtmsd,'("( F12.3, ",I0,"ES15.7 )")') 5*(nmoltype+1)
    if (bCalcTransfer) then
        write(strfmtconduct,'("( F12.3, ",I0,"ES15.7 )")') (nmoltype)*(nmoltype)+1
    else
        write(strfmtconduct,'("( F12.3, ES15.7 )")') 
    endif
    ! write (0,1) for auto correlation functions
    write(20,strfmt) 0.d0,rotACFt0
    write(21,strfmt) 0.d0,rotACFt0
    do idt=1,ndframe
     ! icnt = nACFTime(idt)
      icnt = nDiffTime(idt)
      if(icnt .eq. 0) then
        continue
      endif
!      if(xqAtomsCFTime(idt) .ne. 0) then
!        write(18,'(2E15.7)') dt0*idt,xqAtomsCFTime(idt)/icnt
!      endif
      if(msdTime(idt,5,nmoltype+1) .ne. 0) then
        write(18,strfmtmsd) dt0*idt,msdTime(idt,:,:)/icnt
      endif
      if(xqcomCFTime(idt,idxlast) .ne. 0) then
        write(19,strfmtconduct) dt0*idt,xqcomCFTime(idt,:)/icnt
      endif

      if(rotACFTime(idt,nmoltype+1) .ne. 0) then
        write(20,strfmt) dt0*idt,rotACFTime(idt,:)/icnt
        rotACFTimeP2(idt,:) = 1.5d0*(rotACFTimeP2(idt,:)/icnt)-0.5d0
        write(21,strfmt) dt0*idt,rotACFTimeP2(idt,:)
      endif
    enddo


    write(6,*) 'finished writing output files'
    close(6)
    close(18)
    close(19)
    close(20)
    close(21)

end program calc_xqCF_rotACF
