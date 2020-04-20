!**************************************************
! these are global variables used by many subroutines throughout
! the program
!**************************************************

MODULE variables

  integer :: nmolcat,nmolani,nxtc,nmolsys,nxtcfile,nmoltype
  integer,DIMENSION(:,:,:),ALLOCATABLE:: idxClustIons
  integer,DIMENSION(:,:),ALLOCATABLE:: idxPairIon

  real*8  :: temperature
  real*8,DIMENSION(:,:),ALLOCATABLE:: distPairIon
!  real*8,dimension(:),allocatable :: timestamp
  real*8,dimension(:),allocatable :: corr_IP,corr_IP_cat,corr_IP_ani
  real*8,dimension(:),allocatable :: corr_CP,corr_CP_cat,corr_CP_ani
  real*8,dimension(3) :: cos_vec
!  real*8 :: dt,t0
  character(len=256)         :: strInFile, strOutFile, strTopFile
  character(len=256)         :: strCOSFile, strCOSSqFile
  character(len=256)         :: strIPFile, strCPFile,strConductFile,strRotACFFile,strRotACFP2File
  character(len=256),allocatable,dimension(:) :: strXtcfiles

END MODULE variables
