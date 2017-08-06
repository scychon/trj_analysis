!****************************************************************
! 
! this module reads topology file and save it to the type topfile
!
!****************************************************************

MODUlE topol


    implicit none
    private
 
!   *** topfile type
!   box     - triclinic pbc box of the configuration
!   NATOMS  - number of atoms in the configuration.
!   pos     - positions read in (3,NATOMS)
!   prec    - precision of the coordinates read in
!   STEP    - step number of configuration.
!   STAT    - status of operation. 0 = good
!   time    - time of the configuration
!   xd      - pointer from libxdrfile.
 
!   Should always call init first. Then call read in a loop and do your
!   calculations. After the loops call close.
 
    type, public :: molecule
      Character(5) :: moltypename
      integer      :: numatom, nummol
      real*8       :: molarmass, molcharge
      real*8, allocatable :: atommass(:)
      real*8, allocatable :: atomcharge(:)
    end type
 
    type, public :: topfile
      type(molecule), allocatable :: moltype(:)
!      type(molecule), pointer, allocatable :: sysmol(:)
      integer      :: numsysmol, nummoltype, nsize_mol, nsize_atom
    contains
      procedure :: init => init_top
    end type

 
 
 
contains
 
    ! our wrappers for the topfile class
    subroutine init_top(top,filename_in)
 
        implicit none
        class(topfile), intent(inout) :: top
        character (len=*), intent(in) :: filename_in
        character (len=206) :: filename
        character(50):: line
        character(5):: atomname,molname
        logical :: ex
        integer :: i,j,k,nmoltype,nmol,natom,nmolsys
        real*8  :: mass, molmass, charge
 
        inquire(file=trim(filename_in),exist=ex)
 
        if (ex .eqv. .false.) then
            write(0,*)
            write(0,'(a)') " Error: "//trim(filename_in)//" does not exist."
            write(0,*)
            stop
        end if
 
        ! Set the file name to be read in for C.
        filename = trim(filename_in)
        OPEN(unit=7,file=filename,status='old')
 
        ! Get number of moltypes
        read(7, '(A)') line
        read(7,*) nmoltype
        top % nummoltype = nmoltype

        allocate(top % moltype(nmoltype))

        ! Read atomic masses for each molecule type
        do i=1, nmoltype
          read(7,'(A)') molname
          read(7,*) natom
          top % moltype(i) % moltypename = molname
          top % moltype(i) % molarmass = 0
          top % moltype(i) % molcharge = 0
          top % moltype(i) % numatom = natom
          allocate(top % moltype(i) % atommass(natom))
          allocate(top % moltype(i) % atomcharge(natom))

          do j=1, natom
            read(7,*) atomname, mass, charge
            top % moltype(i) % atommass(j) = mass
            top % moltype(i) % atomcharge(j) = charge
            top % moltype(i) % molarmass = top % moltype(i) % molarmass + mass
            top % moltype(i) % molcharge = top % moltype(i) % molcharge + charge
          enddo
        enddo

        ! Read total # of molecules in system
        read(7, '(A)') line
        read(7, *) nmol
        top % numsysmol = nmol
!        allocate(top % sysmol(nmol))
        nmolsys = nmol

        ! Read # of molecules for each moltype
        k = 0
        do i=1, nmoltype
          read(7,*) molname, nmol
          top % moltype(i) % nummol = nmol
          write(*,*) 'molcharge', top % moltype(i) % molcharge
          write(*,*) 'molarmass', top % moltype(i) % molarmass
          do j=1, nmol
            k = k+1
!            top % sysmol(k) => top % moltype(i)
          enddo
        enddo

        if(k .ne. nmolsys) then
            write(0,*)
            write(0,'(a)') " Error: number of molecules of each types does not sum up to the total number of molecules in system."
            write(0,*)
            stop
        end if

        close(7)

    end subroutine init_top

!    recursive subroutine read_itp(top, itpfile)
!        implicit none
!        class(topfile), intent(inout) :: top
!        character (len=*), intent(in) :: itpfile
!        character (len=206) :: filename
!        character(50):: line
!        character(5):: atomname,molname
!        logical :: ex
!        integer :: i,j,k,nmoltype,nmol,natom,nmolsys
!        integer :: readstat
!        real*8  :: mass, molmass
!        real*8, allocatable :: temp(:)
! 
!        inquire(file=trim(itpfile),exist=ex)
! 
!        if (ex .eqv. .false.) then
!            write(0,*)
!            write(0,'(a)') " Error: "//trim(itpfile)//" does not exist."
!            write(0,*)
!            stop
!        end if
! 
!        ! Set the file name to be read in for C.
!        filename = trim(itpfile)
!        OPEN(unit=7,file=filename,status='old')
!
!        readstat = 0 
!
!        ! Read until end of file
!        do while (readstat .ge. 0)
!          read(15, '(A)',IOSTAT=readstat) line
!
!          !remove comments
!
!          !separate line into variables
!
!          !include statements
!
!          !define statements
!
!          !section titles
!
!
!          !read each section parameters
!
!
!          !
!          if (readstat .ge. 0) then
!            i=i+1
!            if(readstat .eq. 0) then
!              j=j+1
!              if(size(re_to_e_list) .lt. j) then
!                allocate(temp(size(re_to_e_list)+nsize))
!                temp=0
!                temp(:size(re_to_e_list))=re_to_e_list
!                deallocate(re_to_e_list)
!                allocate(re_to_e_list(size(temp)))
!                re_to_e_list=temp
!                deallocate(temp)
!              endif
!              re_to_e_list(j)=re_to_e
!              if(re_to_e .gt. rmax) rmax=re_to_e
!            else
!              write(*,*) 'something is wrong with %ith data', i
!            endif
!          else
!            write(*,*) 'reading file finished'
!          endif
!        enddo
!
!        ! Get number of moltypes
!        read(7, '(A)') line
!    end subroutine

end module topol
