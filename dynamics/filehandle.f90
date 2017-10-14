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
        bMSDVec=.false.
        bConduct=.false.
        bMSD=.false.
        bRotACF=.false.
        bRotACFP2=.false.
        bVelXtc=.false.
        bSingleChain=.true.
        bPolyStat=.false.
        bSq=.false.
        bSqTime=.false.
        bTemp=.false.
        bXtcFiles=.false.

        !loop across options
        do iarg=1,narg
          str = aStrArgs(iarg)
          select case(adjustl(str))
            case("--help","-h")
               write(*,*)"This is program polystat : Version 0.1"
               write(*,*) 'Options : '
          
            !First known args
            case("-f")
               bLookForInp=.TRUE. !change logical value
            case("-o")
               bLookForOut=.TRUE.
            case("-msd")
               bMSD=.TRUE.
            case("-rot")
               bRotACF=.TRUE.
            case("-poly")
               bPolyStat=.TRUE.
            case("-sq")
               bSq=.TRUE.
            case("-sqt")
               bSqTime=.TRUE.
            case("-cond")
               bConduct=.TRUE.
            case("-vel")
               bVelXtc=.TRUE.
            case("-all")
               bMSD=.TRUE.
               bRotACF=.TRUE.
               bPolyStat=.TRUE.
               bConduct=.TRUE.
               bVelXtc=.TRUE.
               bSq=.TRUE.
               bSqTime=.TRUE.
          
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
          write(*,*) 'usage : polyanal -f param.dat [options]'
          stop
        endif
        
        inquire(file=strInFile,exist=bFileExist)!check if it exist
        if(.not.bFileExist)then
          call printerr('file ' // strInFile // ' not found')
        elseif( (.not.bMSD) .and. (.not.bRotACF) .and. (.not.bConduct) .and. (.not.bVelXtc) .and. (.not.bPolyStat) ) then
          call printerr('No analysis option is selected')
        endif
        
        if(strOutFile .eq. "") then
          write(*,*) 'usage : conductivity -f param.dat -o outfile.dat'
          write(*,*) 'output file has not set'
          write(*,*) 'will use default outfile name conductivity_out.dat'
          strOutFile = 'conductivity_out.dat'
        endif

    end subroutine getargs

    ! subroutine to print error message
    subroutine printerr(strError)
        use variables
        implicit none
        character(len=*), intent(in)  :: strError
        write(*,*) 'Error : ',strError
        write(*,*) 'usage : polyanal -f param.dat [options]'
        write(*,*) 'To see the available options, type polyanal -h'
        stop
    end subroutine printerr

    ! subroutine to read the parameter file
    subroutine readparam
        use variables
        use fvector
        implicit none
        integer :: ifile,iline, ISTAT, idx
        character(len=256)         :: str, line, strhead

        strTopFile = ""
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
            cycle   ! no info in the line
          endif

          if (line(:1) .eq. '[') then
            strhead = trim(adjustl(line(2:(index(line,']')-1))))
          else
            write(*,*) 'Unkown parameter on ',iline, 'th line : '
            write(*,*) line
!            write(*,*) 'Error in input file ',strInFile
!            write(*,*) iline, 'th line does not contain proper data !'
!            stop
          endif

          select case (strhead)
            case ("xtcfiles")
              read(7,*) nxtcfile
              write(*,*) nxtcfile
              allocate(strXtcfiles(nxtcfile))
              if (nxtcfile .ge. 0) then
                bXtcFiles = .true.
              endif
              do ifile=1,nxtcfile
                read(7,'(A)') strXtcfiles(ifile)
                write(*,*) trim(strXtcfiles(ifile))
              enddo
            case ("velfiles")
              if (.not.bVelXtc) cycle
              read(7,*) nvelxtcfile
              write(*,*) nvelxtcfile
              allocate(strVelXtcfiles(nvelxtcfile))
              if (nvelxtcfile .ge. 0) then
                bVelXtc = .true.
              endif
              do ifile=1,nvelxtcfile
                read(7,'(A)') strVelXtcfiles(ifile)
                write(*,*) trim(strVelXtcfiles(ifile))
              enddo
            case ('topfile')
              read(7, '(A)') strTopFile
              write(*,*) strTopFile
            case ('polyfile')
              !if (.not.bPolyStat) cycle
              read(7, '(A)') strPolyStatFile
              write(*,*) strPolyStatFile
            case ('msdfile')
              !if (.not.bMSD) cycle
              read(7, '(A)') strMSDFile
              write(*,*) strMSDFile
            case ('msdvec')
              !if (.not.bMSD) cycle
              read(7, *) msd_vec
              write(*,*) msd_vec
              msd_vec = msd_vec/norm(msd_vec)
              bMSDVec = .true.
            case ('conductivity output file')
              !if (.not.bConduct) cycle
              read(7, '(A)') strConductFile
              write(*,*) strConductFile
            case ('rotational ACF file')
              !if (.not.bRotACF) cycle
              read(7, '(A)') strRotACFFile
              write(*,*) strRotACFFile
            case ('rotational ACF P2 file')
              !if (.not.bRotACF) cycle
              read(7, '(A)') strRotACFP2File
              write(*,*) strRotACFP2File
              bRotACFP2 = .true.
            case ('temp')
              read(7, *) temperature
              write(*,*) temperature
              bTemp = .true.
            case ('skip')
              read(7, *) nskip
              write(*,*) 'calculate properties from every ',nskip,'th snapshots only'
            case ('gr_skip')
              read(7, *) nskipgr
              write(*,*) 'calculate g(r) from every ',nskipgr,'th snapshots only'
            case ('gr_bins')
              read(7, *) nbins
              write(*,*) 'number of bins for histogra of g(r) or S(q) : ',nbins
            case ('grfile')
              read(7, *) strGrFile
              write(*,*) strGrFile
            case ('Sqfile')
              read(7, *) strSqFile
              write(*,*) strSqFile
            case ('Sq_rmax')
              read(7, *) sqRmax
              write(*,*) 'calculate g(r) from every ',nskipgr,'th snapshots only'
            case ('AFF_file_dir')
              read(7, *) strAFFDir
              write(*,*) strAFFDir
            case ('AFF_file_ext')
              read(7, *) strAFFExt
              write(*,*) strAFFExt
          end select
          read(7, '(A)',IOSTAT=ISTAT) line   ! search for header line
        enddo
        close(7)

        if(strTopFile .eq. "") then
          write(*,*) 'topology file is not set in your parameter file', strInFile
          write(*,*) 'will use default topfile name topol.dat'
          strTopFile = 'topol.dat'
        endif

        if (.not. bTemp) then
          write(*,*) 'Temperature is not set in your parameter file', strInFile
          write(*,*) 'For help, try polystat_omp -h'
          stop
        elseif (.not. bXtcFiles) then
          write(*,*) 'Input trajectory files are not set in your parameter file', strInFile
          write(*,*) 'For help, try polystat_omp -h'
          stop
        endif

    end subroutine readparam

    subroutine readaff(arrAtomTypes)
        implicit none
        character*4, intent(in), dimension(:)  :: arrAtomTypes
        character(len=256)  :: strFile
        character(len=4)  :: strPreFix
        LOGICAL :: file_exists
        integer :: itype, nAtomType, ISTAT,i

        strPreFix = "AFF_"

        nAtomType = size(arrAtomTypes)
        allocate(atomtype_form(nAtomType, MAX_Q_FORM))
        do itype=1,nAtomType
            strFile = trim(strAFFDir)//'/'//strPreFix//trim(arrAtomTypes(itype))//'.'//trim(strAFFExt)
            OPEN(unit=7,file=strFile,status='old')
            i=0
            do while (ISTAT .eq. 0)
                i = i+1
                if (i .gt. MAX_Q_FORM) then
                    exit
                endif
                read(7, *,IOSTAT=ISTAT) q_grid_form(i), atomtype_form(itype,i)
            enddo
            close(7)
        enddo

        ! convert from nm^-1, to angstrom^-1
        !q_grid_form(:) = q_grid_form(:) / 10d0

    end subroutine

    subroutine calcQr(boxlen)
        implicit none
        integer :: i, j
        real*8 :: qrval, boxlen, dr

        allocate(sinqr_over_qr(max_q_form,nbins))        
        boxlen = boxlen/2.d1    ! calc half the box size in Angstrom unit
        dr = boxlen/nbins
        do i=1, nbins
            do j =1, MAX_Q_FORM
                qrval = q_grid_form(j)*dr*i
                if (qrval < 1e-10) then
                    sinqr_over_qr(j,i) = 1
                else
                    sinqr_over_qr(j,i) = dsin(qrval)/qrval
                endif
            enddo
        enddo
    end subroutine

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
        call system('mv ' // trim(strFile) // ' ' // trim(strFileBak))
        write(*,*) 'file ' // trim(strFile) // ' exists and is backed up to file ' // trim(strFileBak)
      endif
    end subroutine backupfile

end module filehandle

