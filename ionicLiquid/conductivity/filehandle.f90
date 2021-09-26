!****************************************************************
! 
! this module contains utility functions for vector calculation
!
!****************************************************************

MODUlE filehandle
    use variables
    implicit none
 
contains


    pure function getLowerCase(string) result(output)
        ! convert string to lower-case
        character(*), intent(in)    :: string
        integer, parameter          :: DUC = ichar('A') - ichar('a')
        character(len(string))      :: output
        character                   :: ch
        integer                     :: i
        do i = 1,len(string)
            ch = string(i:i)
            if (ch>='A' .and. ch<='Z') ch = char(ichar(ch)-DUC)
            output(i:i) = ch
        end do
    end function getLowerCase

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

    end subroutine getargs

    ! subroutine to read the parameter file
    subroutine readparam
        use variables
        use fvector
        implicit none
        integer :: ifile,iline, ISTAT, idx
        character(len=256)         :: str, line, strhead, strCalcTf

        strTopFile = ""
        bCalcTransfer = .false.
        write(*,*) 'start reading paramdat'
        OPEN(unit=7,file=strInFile,status='old')
        iline=0
        read(7, '(A)',IOSTAT=ISTAT) line   ! search for header line
        do while (ISTAT .eq. 0)
          iline = iline + 1
          line = adjustl(line)
          write(*,*) trim(line)
          ! check whether it has only comments
          idx = min(index(line,'!'),index(line,'#'))
          if (idx .gt. 0) then
            line = line(:idx-1)
          endif
          if (len(trim(line)) .eq. 0) then
            continue   ! no info in the line
          endif

          if (line(:1) .eq. '[') then
            strhead = trim(adjustl(line(2:(index(line,']')-1))))
          else
            write(*,*) 'Error in input file ',strInFile
            write(*,*) iline, 'th line does not contain proper data !'
            stop
          endif

          select case (strhead)
            case ("xtcfiles", "dcdfiles", "trajfiles")
              read(7,*) nxtcfile
              write(*,*) nxtcfile
              allocate(strXtcfiles(nxtcfile))
              do ifile=1,nxtcfile
                read(7,'(A)') strXtcfiles(ifile)
                write(*,*) trim(strXtcfiles(ifile))
              enddo
            case ('topfile')
              read(7, '(A)') strTopFile
              write(*,*) strTopFile
            case ('msdfile')
              read(7, '(A)') strMSDFile
              write(*,*) strMSDFile
            case ('conductivity output file')
              read(7, '(A)') strConductFile
              write(*,*) strConductFile
            case ('rotational ACF file')
              read(7, '(A)') strRotACFFile
              write(*,*) strRotACFFile
            case ('rotational ACF P2 file')
              read(7, '(A)') strRotACFP2File
              write(*,*) iline, strRotACFP2File
            case ('temp')
              read(7, *) temperature
              write(*,*) iline, temperature
            case ('msd_axis','msd_vec')
              read(7, *) msd_vec
              write(*,*) iline, msd_vec
            case ('calc_transference')
              read(7, *) strCalcTf
              write(*,*) iline, strCalcTf
              if (getLowerCase(strCalcTf) .eq. "true") then
                bCalcTransfer = .true.
              endif
          end select
          read(7, '(A)',IOSTAT=ISTAT) line   ! search for header line
        enddo
        close(7)

        if(strTopFile .eq. "") then
          write(*,*) 'topology file is not set in your parameter file ', strInFile
          write(*,*) 'will use default topfile name conparam_bmimbf4.dat'
          strTopFile = 'param_bmimbf4.dat'
        endif

        msd_vec = msd_vec/norm(msd_vec)
        if(strTopFile .eq. "") then
          write(*,*) 'topology file is not set in your parameter file ', strInFile
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

