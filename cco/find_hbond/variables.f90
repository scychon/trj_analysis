!**************************************************
! these are global variables used by many subroutines throughout
! the program
!**************************************************

MODULE variables
    ! 1. Use the xdr interface
    use xtc, only: xtcfile
    ! 2. Declare a variable of type xtcfile
    type(xtcfile) :: traj

  integer :: nmolcat,nmolani,nxtc,nmolsys,natom
  integer :: nneigh, nAtCavity, nGLU, nPRD, nBNC, nRef, idxCuB
  integer,DIMENSION(3):: idxGLU, idxPRD
  integer,DIMENSION(4):: idxBNC
  integer,dimension(:,:),allocatable :: hbondcnt,tempHbond   !matrix for # of water molecules
  integer,dimension(:),allocatable :: tempWatNeigh, idxWatNeigh   !matrix for idxs of water molecules within the cutoff
  integer,DIMENSION(:,:,:,:),ALLOCATABLE :: idxAtomPositions

  logical,dimension(:),allocatable :: bCheckedWater

  real*8             :: dCut, dHbond, angHbond, dCuB
  real*8,dimension(3) :: xgluc, xprdc, xcub          ! position of center of glu oxigens, prd oxigens and cub
  real*8,dimension(:,:),allocatable :: xglu, xprd,xbnc          ! position of center of glu oxigens and prd oxigens
  real*8,dimension(:),allocatable :: timestamp
  real*8,DIMENSION(:,:,:),ALLOCATABLE :: distMat,tempDistMat
  real*8 :: dt,t0

END MODULE variables
