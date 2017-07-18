!///////////////////////////////////////////////////////////////////////
! Initilaize / terminate work areas
module FPinit
  implicit none
contains
!-----------------------------------------------------------------------
subroutine initialize()
  use FPglobal
  use FPvars
  use FPgrid
  implicit none
  integer :: ierror
!
! Setting the lengths of the work areas
!
  NCPTS = KORD*NINT-NCC*(NINT-1)
  ML = NPDE*(KORD+IQUAD-1)-1
  LWORK = KORD+4*NPDE+9*NPDE**2+ &
       NCPTS*(3*KORD+2)+NPDE*NCPTS*(3*ML+MAXDER+7)
  LIWORK = NCPTS*(NPDE+1)
  NBKPT = NINT+1
!
! Assign synonyms
!
  nM = NPDE
  nX = NCPTS
!
  allocate(XBKPT(NBKPT), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate XBKPT(",NBKPT,")"
     stop
  endif
  allocate(WORK(LWORK), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate XBKPT(",NBKPT,")"
     stop
  endif
  allocate(IWORK(LIWORK), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate IWORK(",LIWORK,")"
     stop
  endif
!
  allocate(Xg(nX), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate Xg(",nX,")"
     stop
  endif
  allocate(Ug(nM,nX), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate Ug(",nM,nX,")"
     stop
  endif
  allocate(UXg(nM,nX), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate UXg(",nM,nX,")"
     stop
  endif
  allocate(UXXg(nM,nX), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate UXXg(",nM,nX,")"
     stop
  endif
  allocate(UTg(nM,nX), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate UTg(",nM,nX,")"
     stop
  endif
  allocate(Qg(nM,nX), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate Qg(",nM,nX,")"
     stop
  endif
!
  allocate(x32(nX), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate x32(",nX,")"
     stop
  endif
  allocate(Xm(nM,nX), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate Xm(",nM,nX,")"
     stop
  endif
  allocate(Ym(nM,nX), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate Ym(",nM,nX,")"
     stop
  endif
  return 
end subroutine initialize
!-----------------------------------------------------------------------
subroutine terminate()
  use FPglobal
  use FPvars
  use FPgrid
  implicit none
!
  deallocate(XBKPT)
  deallocate(WORK)
  deallocate(IWORK)
!
  deallocate(Xg)
  deallocate(Ug)
  deallocate(UXg)
  deallocate(UXXg)
  deallocate(UTg)
  deallocate(Qg)
!
  deallocate(x32)
  deallocate(Xm)
  deallocate(Ym)
  return 
end subroutine terminate
end module FPinit
