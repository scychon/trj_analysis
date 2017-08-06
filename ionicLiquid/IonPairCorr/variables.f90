!**************************************************
! these are global variables used by many subroutines throughout
! the program
!**************************************************

MODULE variables

  integer :: nmolcat,nmolani,nxtc,nmolsys
  integer,DIMENSION(:,:,:),ALLOCATABLE:: idxClustIons
  integer,DIMENSION(:,:),ALLOCATABLE:: idxPairIon
  REAL*8,DIMENSION(:,:),ALLOCATABLE:: distPairIon
  real*8,dimension(:),allocatable :: timestamp
  real*8,dimension(:),allocatable :: corr_IP,corr_IP_cat,corr_IP_ani
  real*8,dimension(:),allocatable :: corr_CP,corr_CP_cat,corr_CP_ani
  real*8 :: dt,t0

END MODULE variables
