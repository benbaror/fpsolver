!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Global parameters used by PDECOL 
module FPglobal
  implicit none 
! PDECOL input constants variables
! Note: the collocation grid has NCPTS = KORD*NINT - NCC*(NINT-1) points
!
  integer :: NINT,NPDE,KORD,NCC
!
! Work area size parameters
!
  integer :: IQUAD = 0 ! IF NOT (KORD=3 and null boundary)
!
! ATTENTION:
! PDECOL uses by default the Gauss-Legendre colloacation method
! for NCC = 2, which gives nonsense results (bug? wrong use?). 
! This option has been turned off by setting the parameter NOGAUS=1 
!
  integer, parameter :: MAXDER = 5
  integer, parameter :: NOGAUS = 1
  integer :: NCPTS
  integer :: ML
  integer :: LWORK
  integer :: LIWORK
!
!     Synonyms: The variable notation of the PDE solver expresses 
!     the numerical scheme, whereas that of the driver expresses 
!     the physical problem. The physics are clearer when descriptive 
!     variable names are used. Therefore, the  following synonyms 
!     are adopted in the driver:
!
!     Physics: nX   : No. of dimensionless energy grid points
!     =
!     PDE    : NCPTS: No. of collocation grid points
!
!     Physics: nM   : No. of mass groups
!     =
!     PDE    : NPDE : No. of PDEs
!
!     Physics: g    : DF (one for each mass group)
!     =
!     PDE    : U    : Function(s) of x to be solved
!
!     Physics: gX   : dg/dx
!     =
!     PDE    : UX   : dU/dx
!
!     Physics: gXX  : d^2g/dx^2
!     =
!     PDE    : UXX  : d^2U/dx^2
!
!     Physics: xU  : The zero-energy boundary (x=0)
!     =
!     PDE    : xL   : The left hand side boundary
!
!     Physics: xD  : The stellar destruction energy boundary
!     =
!     PDE    : xR   : The right hand side boundary
!
  integer :: nX, nM
!
end module FPglobal
