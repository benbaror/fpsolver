!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Global parameters used by PDECOL 
module FPvars
!prm!  use FPglobal
  implicit none 
!
!     Version string
!
  character(len=10), parameter :: version =  'V080627g95'
!
! Work areas
!
  integer :: NBKPT
!prm!  integer, parameter :: NBKPT = NINT+1
!
  real, dimension(:), allocatable, save :: XBKPT
  real, dimension(:), allocatable, save :: WORK
  integer, dimension(:), allocatable, save :: IWORK

!
! Additional PDECOL variables
!
  integer :: MF
  real :: DT, EPS
!
! Description of physical system: mass, number ratios, boundary conds
!
  integer, parameter :: mM = 64 ! needs to be predefined because 
                                ! masses and boundary conds
                                ! appear in namelist
  real, dimension(mM), save :: m, Cm, Dm, Sm
  REAL ::  xU,xD,x1
  REAL :: ARR, ANR, Aeps, AGR, xmin, xmax, eta_L
  real :: rh2rg                ! ratio of rh=GM/sig(*)^2 to rg=GM/c^2
                               ! used for calculating LC term
  integer :: iMtyp             ! Index of typical star
!new(
!
! Read initial conditions from a file (if 1st char != " ")
!
  character(len=100), save :: initf = ' '
  integer, save :: x_idx, g_idx(mM)
!new)
!
! Debug flags
!
  logical :: DEBUG
  integer :: DBGLVL
!
! A flagging file: if it exists, the program will execute a clean
! exit and write the results as if terminating normally.
!
  character(len=9), parameter :: stopfl = 'stop.run'
!
! Vectors of constants for the calculation of the currents Q that
! need to be calculated only once (by subrouting Qginit) *after* the 
! collocation grid has been initialized by PDECOL.
!
! x32=x^(3/2), Xm=m/x^(3/2), Ym=m^2*Cm/Dm/x^(3/2)
  real, dimension(:), allocatable :: x32
  real, dimension(:,:), allocatable :: Xm
  real, dimension(:,:), allocatable :: Ym
! Flag to indicate whether Qgcons variables need to be initialized
  LOGICAL :: QCinit = .TRUE.
  REAL, PARAMETER :: pi = 4*ATAN(1.0)
!
end module FPvars
