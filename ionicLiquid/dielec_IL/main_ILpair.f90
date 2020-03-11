!  XDR Fortran Interface xtc Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!
 
program calc_dielect
 
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
  integer            :: nframe


  logical::lookForInp=.FALSE.
  logical::lookForOut=.FALSE.
  logical::fileExist

  real*8             :: delr=0.02d0,rmax=3.0d0,boxsize=0.3d0      !
  real*8             :: mass, molmass, dr, rclustcat, rclustani,charge,molcharge
  real*8             :: murotsq, mutranssq, murottrans, volume,temperature
  real*8             :: rotsq,transsq,rottrans,rotval,transval
  real*8             :: debye, eps0, kb, enm, multiple, qelec,mult2
  real*8             :: volavg, muavg
  real*8,dimension(3)    :: xpos,sysbox,compos,dist,box,mutot
  real*8,dimension(3)    :: murot, mutrans
  real*8,dimension(300)  :: rbin
  real*8, allocatable,dimension(:,:,:) :: comtraj,comtrajcat,comtrajani,mutraj,mutrajcat,mutrajani,temp   ! trajectory of com of each molecules (frame idx, mol idx, xyz)
  real*8, allocatable,dimension(:,:) :: boxtraj,tempbox      !matrix for # of water molecules along z direction
  real*8, allocatable,dimension(:,:) :: mutrajrot, mutrajtrans      !matrix for murot,mutran
  real*8, allocatable,dimension(:) :: temptime     !matrix for # of water molecules along z direction
  character(len=256)         :: strout,strtitle
  character(len=256)         :: strInFile, strOutFile, strTopFile
  character(len=256)         :: strDielecFile, strIPFile, strCPFile
  character(len=256)         :: str, line
  character(len=256),allocatable,dimension(:) :: strXtcfiles
! ----------------------------------------------------------------

  murotsq = 0
  mutranssq = 0
  murottrans = 0
  muavg = 0
  nframe = 0
  volavg = 0
  volume = 0

  debye = 3.33564E-30
  eps0 = 8.85418781762E-12
  kb = 1.380648813E-23
  enm = 0.020819434
  qelec = 1.60217656535E-19

  multiple = debye**2/(enm**2 * eps0 * kb)
  mult2 = qelec**2 * 1e-18/(eps0*kb)

  write(*,*) multiple, mult2

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
  write(*,*) 'usage : dielectric -f param.dat -o outfile.dat -n index.ndx'
  stop
endif

inquire(file=strInFile,exist=fileExist)!check if it exist
if(.not.fileExist)then
  write(*,*)'file ',strInFile,' not found'
  write(*,*) 'usage : dielectric -f param.dat -o outfile.dat -n index.ndx'
  stop
endif

if(strOutFile .eq. "") then
  write(*,*) 'usage : dielectric -f param.dat -o outfile.dat -n index.ndx'
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
    read(7, *) strDielecFile
    read(7, '(A)') line
    read(7, *) temperature
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
    allocate(mutraj(nsize,nmolsys,3))
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

    OPEN (UNIT=6,FORM='FORMATTED')
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
        volume = traj % box(1,1) * traj % box(2,2) * traj % box(3,3)
        volavg = volavg + volume
!        write(*,*) 'volume', volume, traj % box(1,1), traj % box(2,2), traj % box(3,3)

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

          allocate(temp(isize+nsize,nmolsys,3))
          temp=0
          temp(:isize,:,:)=mutraj
          deallocate(mutraj)
          allocate(mutraj(isize+nsize,nmolsys,3))
          mutraj=temp
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
            mutot = 0

            do k=1, sys % moltype(i) % numatom
              idxatom = idxatom +1
              mass = sys % moltype(i) % atommass(k)
              charge = sys % moltype(i) % atomcharge(k)
              xpos(:) = traj % pos(:,idxatom)
              compos(:) = compos(:) + xpos(:)*mass
              mutot(:) = mutot(:) + xpos(:)*charge
            enddo

            molmass = sys % moltype(i) % molarmass
            compos(:) = compos(:)/molmass
            comtraj(ixtc,idxmol,:) = compos(:)
            mutraj(ixtc,idxmol,:) = mutot(:)
!            write(*,*) comtraj(ixtc,idxmol,:)
          enddo
        enddo 

        write(6,100,advance='no') achar(13), ixtc,'th frame has finished  '
 100    FORMAT(A,I5,A)

        call traj % read
 
    end do

    ! 5. Close the file
    call traj % close
    end do

    nxtc = ixtc
    nframe = nframe + nxtc
    volavg = volavg * 1e-27 / nxtc
    multiple = multiple / (temperature * volavg * 3)
    write(*,*) 'temp mult vol', multiple, temperature, volavg

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

      allocate(temp(nxtc,nmolsys,3))
      temp=0
      temp(:,:,:)=mutraj(:nxtc,:,:)
      deallocate(mutraj)
      allocate(mutraj(nxtc,nmolsys,3))
      mutraj=temp
      deallocate(temp)
    endif

    write(6,*) 'copying position data into each ion matrices' 
    
    allocate(comtrajcat(nxtc,nmolcat,3))
    allocate(comtrajani(nxtc,nmolani,3))
    allocate(mutrajcat(nxtc,nmolcat,3))
    allocate(mutrajani(nxtc,nmolani,3))
    allocate(mutrajrot(nxtc,3))
    allocate(mutrajtrans(nxtc,3))
    allocate(idxPairIon(nxtc,nmolsys))
    allocate(distPairIon(nxtc,nmolsys))
    allocate(idxClustIons(nxtc,nmolsys,10))
    comtrajcat(:,:,:) = comtraj(:,:nmolcat,:)
    comtrajani(:,:,:) = comtraj(:,(nmolcat+1):,:)
    mutrajcat(:,:,:) = mutraj(:,:nmolcat,:)
    mutrajani(:,:,:) = mutraj(:,(nmolcat+1):,:)
    mutrajrot = 0
    mutrajtrans = 0

    write(6,*) 'start generating neighborlist data'
    idx = 0
    write(6,*) 'first frame start'
    do ixtc=1, nxtc
      box(:) = boxtraj(ixtc,:)
      dr = norm(box*0.5d0)
      distPairIon(ixtc,:) = dr

      idxmol = 0
      do i=1, sys % nummoltype
        do j=1, sys % moltype(i) % nummol
          idxmol = idxmol + 1
          molcharge = sys % moltype(i) % molcharge
          mutrans = comtraj(ixtc,idxmol,:)*molcharge
          mutrajtrans(ixtc,:) = mutrajtrans(ixtc,:) + mutrans
          mutrajrot(ixtc,:) = mutrajrot(ixtc,:) + mutraj(ixtc,idxmol,:) - mutrans
        if (idxmol .le. 200) then
          muavg = muavg + norm(mutraj(ixtc,idxmol,:) - mutrans)
        endif
        enddo
      enddo

      murotsq =  murotsq + dot_product(mutrajrot(ixtc,:),mutrajrot(ixtc,:))
      mutranssq =  mutranssq + dot_product(mutrajtrans(ixtc,:),mutrajtrans(ixtc,:))
      murottrans =  murottrans + dot_product(mutrajrot(ixtc,:),mutrajtrans(ixtc,:))
      write(6,100,advance='no') achar(13), ixtc,'th frame has finished  ' 
!      write(*,*) ixtc,dot_product(mutrajrot(ixtc,:),mutrajrot(ixtc,:)),dot_product(mutrajtrans(ixtc,:),mutrajtrans(ixtc,:)) 
    enddo

    murotsq = murotsq/nxtc
    mutranssq = mutranssq/nxtc
    murottrans = murottrans/nxtc
    muavg = muavg/(nxtc*nmolcat)
    

!    mutrajrot=mutrajrot*multiple
!    mutrajtrans=mutrajtrans*multiple
    murotsq = murotsq*multiple
    mutranssq = mutranssq*multiple
    murottrans = murottrans*multiple
    muavg = muavg /enm

    open(18,file=strDielecFile)
    rotsq = 0
    transsq = 0
    rottrans = 0
    murot = 0
    mutrans = 0
    do i=1,nxtc
      write(18,'(I5,3E15.7)') i,dot_product(mutrajrot(i,:),mutrajrot(i,:)), &
        dot_product(mutrajtrans(i,:),mutrajtrans(i,:)),dot_product(mutrajrot(i,:),mutrajtrans(i,:))
      rotsq = rotsq + dot_product(mutrajrot(i,:),mutrajrot(i,:))
      transsq = transsq + dot_product(mutrajtrans(i,:),mutrajtrans(i,:))
      rottrans = rottrans + dot_product(mutrajrot(i,:),mutrajtrans(i,:))
      murot(:) = murot(:) + mutrajrot(i,:)
      mutrans(:) = mutrans(:) + mutrajtrans(i,:)
    enddo
    rotsq = rotsq/nxtc
    transsq = transsq/nxtc
    rottrans = rottrans/nxtc
    murot = murot/nxtc
    mutrans = mutrans/nxtc
    write(18,'(A5,3E15.7)') 'tot  ',rotsq,transsq,rottrans
    write(18,'(A5,3E15.7)') 'tot  ',dot_product(murot,murot),dot_product(mutrans,mutrans)
    write(18,'(A15,2f8.4)') 'dielectric  ',murotsq+1, 1+murotsq+murottrans
    write(18,'(A15,2f8.4)') 'dielectric  ',rotsq*multiple+1, 1+rotsq*multiple+rottrans*multiple
    write(18,'(A15,1f8.4)') 'avgmu  ',muavg
    close(18) 

    write(*,*) murotsq, mutranssq, murottrans
    write(*,*) rotsq, transsq, rottrans
    write(*,'(A5,3E15.7)') 'tot  ',dot_product(murot,murot),dot_product(mutrans,mutrans)
    write(*,'(A15,2f8.4)') 'dielectric  ',murotsq+1, 1+murotsq+murottrans
    write(*,'(A15,2f8.4)') 'dielectric  ',rotsq*multiple+1, 1+rotsq*multiple+rottrans*multiple
    write(*,'(A15,1f8.4)') 'avgmu  ',muavg
    write(*,'(A15,1f8.4)') 'multiple  ',multiple



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

end program calc_dielect
