!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FPSolver grid areas (for implementing the differential-integral PDEs)
module FPgrid
  implicit none 
! Xg(nX), Ug(nM,nX), UXg(nM,nX), UXXg(nM,nX), UTg(nM,nX), Qg(nM,nX)
  real, dimension(:), allocatable :: Xg
  real, dimension(:,:), allocatable :: Ug
  real, dimension(:,:), allocatable :: UXg
  real, dimension(:,:), allocatable :: UXXg
  real, dimension(:,:), allocatable :: UTg
  real, dimension(:,:), allocatable :: Qg
  integer :: iXg
end module FPgrid
