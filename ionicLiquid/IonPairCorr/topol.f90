!****************************************************************
! 
! this module reads topology file and save it to the type topfile
!
!****************************************************************

module topol


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
      real*8       :: molarmass
      real*8, allocatable :: atommass(:)
    end type
 
    type, public :: topfile
      type(molecule), allocatable :: moltype(:)
!      type(molecule), pointer, allocatable :: sysmol(:)
      integer      :: numsysmol, nummoltype
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
        real*8  :: mass, molmass
 
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
          top % moltype(i) % numatom = natom
          allocate(top % moltype(i) % atommass(natom))

          do j=1, natom
            read(7,*) atomname, mass
            top % moltype(i) % atommass(j) = mass
            top % moltype(i) % molarmass = top % moltype(i) % molarmass + mass
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

end module topol
