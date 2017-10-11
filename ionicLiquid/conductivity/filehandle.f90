!****************************************************************
! 
! this module contains utility functions for vector calculation
!
!****************************************************************

MODUlE filehandle
    use variables
    implicit none
 
contains

    ! Subroutine to get command arguments and assign input and output file name
    subroutine getargs ( aStrArgs )
        use variables
        implicit none

        character(len=256),intent(in),allocatable,dimension(:) :: aStrArgs
        integer :: narg,iarg
        logical::bLookForInp=.FALSE.
        logical::bLookForOut=.FALSE.
        logical::bFileExist
        character(len=256)         :: str, line

        ! Initialization
        strInFile = ""
        strOutFile = ""
        narg=size(aStrArgs)
        !loop across options
        do iarg=1,narg
          str = aStrArgs(iarg)
          select case(adjustl(str))
            case("--help","-h")
               write(*,*)"This is program TestArg : Version 0.1"
          
            !First known args
            case("-f")
               bLookForInp=.TRUE. !change logical value
            case("-o")
               bLookForOut=.TRUE.
          
            case default
            !Treat the second arg of a serie
              if(bLookForInp)then
                strInFile=adjustl(str) !assign a value to pedfile
                bLookForInp=.FALSE. !put the logical variable to its initial value
              elseif(bLookForOut)then
                strOutFile=adjustl(str)
                inquire(file=strOutFile,exist=bFileExist)
                if(bFileExist)then
                 write(*,*)'file ',strOutFile,' exist'
                endif
                bLookForOut=.FALSE.
              else
                write(*,*)"Option ",adjustl(str),"unknown"
              endif
          end select
        end do

        if(strInFile .eq. "") then
          write(*,*) 'input file has not set'
          write(*,*) 'usage : conductivity -f param.dat -o outfile.dat'
          stop
        endif
        
        inquire(file=strInFile,exist=bFileExist)!check if it exist
        if(.not.bFileExist)then
          write(*,*)'file ',strInFile,' not found'
          write(*,*) 'usage : conductivity -f param.dat -o outfile.dat'
          stop
        endif
        
        if(strOutFile .eq. "") then
          write(*,*) 'usage : conductivity -f param.dat -o outfile.dat'
          write(*,*) 'output file has not set'
          write(*,*) 'will use default outfile name conductivity_out.dat'
          strOutFile = 'conductivity_out.dat'
        endif

    end subroutine getargs

    ! subroutine to read the parameter file
    subroutine readparam
        use variables
        use fvector
        implicit none
        integer :: ifile
        character(len=256)         :: str, line

        strTopFile = ""
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
        write(*,*) line
        read(7, '(A)') strTopFile
        write(*,*) strTopFile
        read(7, '(A)') line
        read(7, '(A)') strMSDFile
        read(7, '(A)') line
        read(7, '(A)') strConductFile
        read(7, '(A)') line
        read(7, '(A)') strRotACFFile
        read(7, '(A)') line
        read(7, '(A)') strRotACFP2File
        read(7, '(A)') line
        read(7, *) temperature
        read(7, '(A)') line
        read(7, *) msd_vec
        close(7)

        msd_vec = msd_vec/norm(msd_vec)
        if(strTopFile .eq. "") then
          write(*,*) 'topology file is not set in your parameter file', strInFile
          write(*,*) 'will use default topfile name conparam_bmimbf4.dat'
          strTopFile = 'param_bmimbf4.dat'
        endif

    end subroutine readparam

    subroutine backupfile ( strFile )
      implicit none
      character(len=256), intent(in)   :: strFile
      character(len=256)               :: strFileIn, strFileBak
      LOGICAL :: file_exists
      integer :: i
      INQUIRE(FILE=strFile, EXIST=file_exists)
      if(file_exists) then
        strFileBak = trim(strFile) // '.bak'
        strFileIn = strFileBak
        i = 0
        INQUIRE(FILE=strFileBak, EXIST=file_exists)
        do while(file_exists)
          write(strFileBak,'(A,I0)') trim(strFileIn) // '.',i
          i = i+1
          INQUIRE(FILE=strFileBak, EXIST=file_exists)
        enddo
        call system('cp ' // trim(strFile) // ' ' // trim(strFileBak))
        write(*,*) 'file ' // trim(strFile) // ' exists and is backed up to file ' // trim(strFileBak)
      endif
    end subroutine backupfile

end module filehandle

