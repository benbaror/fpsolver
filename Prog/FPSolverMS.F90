! Driver for FPSolverMS
!///////////////////////////////////////////////////////////////////////
program FPSolverMS
  use FPglobal
  use FPvars
  use FPinit
  use FPwork
  use FPgrid
  use FPcurrent
  use FPpdecol
  implicit none
!
! PDECOL parameters
!
  integer :: INDEX
  real :: T0, T1
!
! Areas for solution 
!
  integer, parameter :: NDERV = 2, NDERV1 = NDERV+1
  integer :: LSCTCH
  integer :: NSOL
  real, dimension(:), allocatable, save :: XSOL
  real, dimension(:,:,:), allocatable, save :: USOL
  real, dimension(:), allocatable, save :: SCTCH
!
! Work areas
!
  real, dimension(:,:), allocatable :: WORKg1, WORKg2
  real, dimension(:,:), allocatable :: UVAL
!
  real, dimension(:), allocatable :: dQ, sQ, eQ
  real :: crate0
  real, dimension(:,:), allocatable :: Ug0
!
  integer :: i,iT,k,id,im,ioerr,nT1
  real :: d,dT1,crate
  character(len=100) :: fname, xgridf = " "
  CHARACTER(len=100) :: run_title = " "
  CHARACTER(len=100) :: path = "./"
  integer :: ierror,iostat
  real :: uti, utf
  character :: fdate*24
!
  namelist /FPinp/ NINT, NPDE, KORD, NCC, &
     DEBUG,DBGLVL,crate0,nT1,dT1,xD,rh2rg, &
     x1,xgridf,initf,x_idx,g_idx,DT,EPS,MF,m,Cm,Dm,Sm, &
     run_title, ANR, ARR, Aeps, AGR, path
!
  write (*,'("Program FPSolver (",a10,") begins.")') version
!
! Verifying that stop flag is not set 
!
!!  if (access(stopfl,'r') .eq. 0) then
!!     write (*,*) "stop.run file exists. Delete it and re-run."
!!     stop
!!  else 
!!     write (*,'("For a clean stop, create the flagging file ",a8)') stopfl
!!  endif
  write (*,'("For a clean stop, create the flagging file ",a8)') stopfl
  write (*,*)
!
! Read input parameters
!
! Initialize boundary conditions to "null" values
!
  Sm = -1.0
  Cm = -1.0  
  Dm = +1.0
!
  read (*,nml=FPinp)
  if (DEBUG) then
     write (*,nml=FPinp)
     write (*,*)
  end if
  if (NPDE > mM) then
     write (*,*) "ERROR: NPDE = ",NPDE, " > mM = ",mM
     stop
  end if
  OPEN (10,file=TRIM(path)//'/FP.inp',status='unknown',iostat=ioerr)
  if (ioerr /= 0) then
     write (*,*) "Error ",ioerr," opening file FP.inp"
     stop
  end if
  write (10,nml=FPinp)
  write (10,'("! Program FPSolver (",a10,")")') version	
#ifdef __LCTERM__
  write(10,'("! With a loss-cone term in the FP equation")')
#else
  write(10,'("! Without a loss-cone term in the FP equation")')
#endif	
#ifdef do_RR
  WRITE(10,'("! With a RR loss-cone term")')
#else
  write(10,'("! Without a RR loss-cone term")')
#endif
#ifdef pinhole
  WRITE(10,'("! Use a pinhole for the full loss-cone")')
#else
  write(10,'("! BW77  loss-cone term")')
#endif

  write (10,'("! ",a)') run_title
  close (10,iostat=ioerr)
  if (ioerr /= 0) then
     write (*,*) "Error ",ioerr," closing file FP.inp"
     stop
  end if
!
! Record start time
!
  uti = usert()
!
! Set paramters and allocate global work areas
!
  call initialize()
!
! Allocate local work areas
!  
  allocate(WORKg1(nX,3), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate WORKg1(",nX,3,")"
     stop
  endif  
  allocate(WORKg2(nX*nM,5), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate WORKg2(",nX*nM,5,")"
     stop
  endif
  allocate(UVAL(nM,3), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate UVAL(",nM,3,")"
     stop
  endif
  allocate(dQ(nM), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate dQ(",nM,")"
     stop
  endif
  allocate(sQ(nM), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate sQ(",nM,")"
     stop
  endif  
  allocate(eQ(nM), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate eQ(",nM,")"
     stop
  endif
  allocate(Ug0(nM,nX), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate Ug0(",nM,nX,")"
     stop
  endif
!
  LSCTCH = KORD*NDERV1
  NSOL = NBKPT
  allocate(XSOL(NSOL), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate XSOL(",NSOL,")"
     stop
  endif
  allocate(USOL(nM,NSOL,NDERV+1), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate USOL(",nM,NSOL,NDERV+1,")"
     stop
  endif
  allocate(SCTCH(LSCTCH), STAT=ierror)
  if(ierror /= 0) then
     write(*,*) "Error trying to allocate SCTCH(",LSCTCH,")"
     stop
  endif
!
  call prthdr(6,dT1,nT1,crate0,run_title,xgridf)
!
!     Initialize work areas
!
!
! Either Sm and Dm are given, and Cm is calculated
! or Cm and Dm are given, and Sm is calculated
!
  where (Sm < 0.0 .and. Cm < 0.0)
     Sm = -1.0
  elsewhere (Sm < 0.0)
     Sm = Cm/Dm**1.5
  elsewhere 
     Cm = Sm*Dm**1.5
  end where
!
  if (minval(Sm(1:nM)) == -1.0) then
     write (*,'("Illegal boundary condition data")')
     stop
!? else 
!?     write (*,'("m  =",1p,9(1x,e8.2))')  m(1:nM)
!?     write (*,'("Sm =",1p,9(1x,e8.2))') Sm(1:nM)
!?     write (*,'("Cm =",1p,9(1x,e8.2))') Cm(1:nM)
!?     write (*,'("Dm =",1p,9(1x,e8.2))') Dm(1:nM)    
  end if
!
! Locate the typical star in the list (Sm~1)
!
  iMtyp = transfer(minloc(abs(Sm-1.0)),0)
  write (*,'("Typical star is m(",i1,") = ",f6.2)') iMtyp,m(iMtyp)
  write (*,*)
!
  WORK = 0.0
  IWORK = 0
  IWORK(1) = LWORK
  IWORK(2) = LIWORK
!
! Create X-grids for PDECOL work grid and solution grid
!
! If 1st character in xgridf is blank, create automatic X-grid
!
      if (xgridf(1:1) == ' ') then    
!
! Create a logarithmically spaced X-grid between some small x1
! and the destruction xD (x=0 is the first point, included for 
! the boundary conditions).
!
         d = (xD/x1)**(1./(NBKPT-2))
         XBKPT(1) = 0.0
         XBKPT(2:NBKPT) = (/ (x1*d**i, i=0,NBKPT-2) /)
!
         d = (xD/x1)**(1./(NSOL-2))
         XSOL(1) = 0.0
         XSOL(2:NSOL) = (/ (x1*d**i, i=0,NSOL-2) /)
      else
!
!     Read X-grid from file (using same work grid and solution grid)
!     Parameters NINT have to be set to #grid-1 and NSOL to #grid 
!      
         if (NSOL .ne. NBKPT) then
            write (*,*) "Error: Manual x steps option, but NSOL <> NBKPT"
            stop
         endif          
         open (10,file=xgridf,status='old',iostat=ioerr)
         if (ioerr /= 0) then
            write (*,*) "Error ",ioerr," opening file ",xgridf
            stop
         end if
         read (10,*) i          ! First line is number of x-values
         if (i .ne. NBKPT) then
            write (*,*) "Error: Manual x steps option, but Nsteps <> NBKPT"
            stop
         endif
         read (10,*) (XBKPT(i),i=1,NBKPT)
         XSOL = XBKPT
         close (10,iostat=ioerr)
         if (ioerr /= 0) then
            write (*,*) "Error ",ioerr," closing file ",xgridf
            stop
         end if
      endif 
      xU = XBKPT(1)
      xD = XBKPT(NBKPT)
!         
      if (DEBUG) write (11,'(i5,1x,1p,e9.3)') (i,XBKPT(i),i=1,NBKPT)
!
      INDEX = 1
      T0 = 0.0
      Ug0 = 0.0
!
!     WARNING: print formats can accomodate only <10 mass types
!
      fname = TRIM(path)//'/FP.out'
      open (10,file=fname,status='unknown',iostat=ioerr)
      if (ioerr /= 0) then
         write (*,*) "Error ",ioerr," opening file ",fname
         stop
      end if
      call prthdr(10,dT1,nT1,crate0,run_title,xgridf)
      fname = 'FP.log'
      open (11,file=fname,status='unknown',iostat=ioerr)
      if (ioerr /= 0) then
         write (*,*) "Error ",ioerr," opening file ",fname
         stop
      end if
!dbg      write (11,'("#",2(a10,1x),27(a9,i1,1x))')
      write (11,'("#",2(a10,1x),192(a9,i1,1x))')
      call prthdr(11,dT1,nT1,crate0,run_title,xgridf)
      write (11,'("#",3(a10,1x), & 
!dbg           & 9("     <Q",i1,">",1x, &
           & 64("     <Q",i1,">",1x, &
           &  "  Delta Q",i1,1x, &
           &  "    rms Q",i1,1x,:))') &
           "T","user time","cnvrg rate",(k,k,k,k=1,nM)
!
      tloop: do iT = 0,nT1
         T1 = T0+dT1*real(iT)
!
         call PDECOL(T0,T1,DT,XBKPT,EPS,NINT,KORD,NCC,nM,MF, &
              INDEX,WORK,IWORK)
         write (*,'("# T1 = ",1p,e9.3," INDEX = ",i3, &
               &  " Elapsed time = ",e9.3)') &
               T1,INDEX,usert()-uti
         if (INDEX .lt. 0) then
            write (*,'(/,"PDECOL error ",i3,". Program aborted.")') INDEX
            stop 
         end if 
!
!     Calculating Q on the collocation points 
!
         do i = 1,nX
            Xg(i) = WORK(IW3+i-1) ! = XC(i)
!     EVAL is an internal PDECOL routine for calculating U,UX,UXX
!     on the collocation grid, which uses data stored in 
!     WORK and IWORK (see explanations in pdecol.f)
            call EVAL(i,NPDE,WORK(IW10),UVAL,WORK(IW1),IWORK)
!                     i NPDE C          UVAL A         ILEFT
            do k = 1,nM
               Ug(k,i)   = UVAL(k,1)
               UXg(k,i)  = UVAL(k,2)
               UXXg(k,i) = UVAL(k,3)               
            end do 
         end do 
         call Qgrid(Qg,Xg,Ug,UXg,UXXg, &
              WORKg1(1,1),WORKg1(1,2), &
              WORKg2(1,1),WORKg2(1,2),WORKg2(1,3),WORKg2(1,4),WORKg2(1,5))
!
! Estimating the constancy of the current by two measures:
! (1) Delta Q = max_x[Q(x)]-min_x[Q(x)]
! (2) rms Q 
!
         call Qcnvrg(Qg,nM,nX,dQ,eQ,sQ)
         call pcnvrg( 6,eQ,dQ,sQ)
         call pcnvrg(10,eQ,dQ,sQ)
!dbg         write (11,'(1p,30(1x,e10.3))') &
         write (11,'(1p,195(1x,e10.3))') &
              T1,usert()-uti,crate,(eQ(k),dQ(k),sQ(k),k=1,nM)
!
         if (DEBUG) then
            do i = 1,nX
!dbg               write (11,'(i3,1p,37(1x,e12.5:))') &
               write (11,'(i3,1p,257(1x,e12.5:))') &
                    iT,Xg(i), &
                    (Ug(k,i),UXg(k,i),UXXg(k,i),Qg(k,i),k=1,nM)
            end do
         end if
!
! Checking the global convergence of the solution
!
         d = sqrt(sum((Ug-Ug0)**2)/(nM*nX))
         crate = d/dT1
         write ( 6,'("# Mean global convergence rate = ",1p,e8.2)') &
              crate
         write ( 6,'("# Elapsed time = ",1p,e9.3)') usert()-uti
         write (10,'("# Mean global convergence rate = ",1p,e8.2)') &
              crate
         write (10,'("# Elapsed time = ",1p,e9.3)') usert()-uti
         if (crate .le. crate0) then
            write ( 6,'("#",/,"# Global convergence achieved at T1 = ", &
                 & 1p,e8.2," : ",1p,e8.2, &
                 & " = crate < crate0 = ",e8.2)') T1,crate,crate0
            write (10,'("#",/,"# Global convergence achieved at T1 = ", &
                 & 1p,e8.2," : ",1p,e8.2, &
                 & " = crate < crate0 = ",e8.2)') T1,crate,crate0
            exit
         end if
!
! Updating the previous step solution
!
         Ug0 = Ug
!
! Getting the solution on the requested (interpolated) solution grid
!
! CAVEAT!!!: Q inexact when calculated on the NINT breakpoints: 
!            Must be calculated on nX collocation points!
!
         call VALUES(XSOL,USOL,SCTCH,nM,NSOL,NSOL,NDERV,WORK)
!
!dbg         write (10,'("#",2(a4,1x),1(a12,1x),28(a11,i1,1x))') &
         write (10,'("# T1 = ",1p,e9.3," x1 = ",e9.3)') T1,x1
         write (10,'("#",2(a4,1x),1(a12,1x),193(a11,i1,1x))') &
              '#X','#T','X',('U',k,'UX',k,'UXX',k,k=1,nM),'Trel'
         do i = 1,NSOL
!dbg            write (10,'(2(1x,i4),1p,1(1x,e12.5),(28(1x,e12.5)))') &
            write (10,'(2(1x,i4),1p,1(1x,e12.5),(193(1x,e12.5)))') &
                 i,iT,XSOL(i),((USOL(k,i,id),id=1,NDERV+1),k=1,nM), &
                 trel(i,iMtyp)
         end do 
         write (10,'(" ")')
!
! Check if a user-requested premature stop is signaled by the
! existence of the file stopfl (defined as a parameter above)
!
         ! if (access(stopfl,'r') .eq. 0) then
         !    open(12,file=stopfl,status='old',iostat=ioerr)
         !    if (ioerr /= 0) then
         !       write (*,*) "Error ",ioerr," opening file ",stopfl
         !       stop
         !    end if
         !    close(12,status='delete',iostat=ioerr)
         !    if (ioerr /= 0) then
         !       write (*,*) "Error ",ioerr," deleting file ",stopfl
         !       stop
         !    end if
         !    write ( 6,'("# Stop signal received")')
         !    write (10,'("# Stop signal received")')
         !    exit
         ! end if
!new(
         open  (12,file='FP.status',status='unknown', &
              position='append',iostat=ioerr)
         if (ioerr /= 0) then
            write (*,*) "Error ",ioerr," opening file FP.status"
            stop
         end if
         write (12,'("T1 = ",1p,e9.3," global convergence rate = ", &
         & e9.3," elapsed time = ",e9.3," date = ",a24)') & 
              T1,crate,usert()-uti,fdate()
         close (12,iostat=ioerr)
         if (ioerr /= 0) then
            write (*,*) "Error ",ioerr," closing file FP.status"
            stop
         end if
!new)
!       
      end do tloop
!
! Convergence status message
!
      if (crate > crate0) then
         write ( 6,'("#",/,"# Global convergence NOT achieved by T1 = ", &
              & 1p,e8.2," : ",1p,e8.2, &
              & " = crate > crate0 = ",e8.2)') T1,crate,crate0
         write (10,'("#",/,"# Global convergence NOT achieved by T1 = ", &
              & 1p,e8.2," : ",1p,e8.2, &
              & " = crate > crate0 = ",e8.2)') T1,crate,crate0
      end if
!     
! Total elapsed time message
!
      write (6,'("# T1 = ",1p,e8.2)') T1
      write (6,'("# Total elapsed time = ",1p,e9.3," s")') &
           usert()-uti
      write (10,'("# T1 = ",1p,e8.2)') T1
      write (11,'("# T1 = ",1p,e8.2)') T1
!
!
      close (11,iostat=ioerr)
      if (ioerr /= 0) then
         write (*,*) "Error ",ioerr," closing file FP.log"
         stop
      end if
      write (10,'("# Total elapsed time = ",1p,e9.3," s")') &
           usert()-uti
      close (10,iostat=ioerr)
      if (ioerr /= 0) then
         write (*,*) "Error ",ioerr," closing file FP.out"
         stop
      end if
      write (*,*)
      write (*,'("# Solutions written on ",a50)') fname
!
!     Writing the last (hopefully steady state) solution
!
      fname = TRIM(path)//'/FP_last.out'
      open (10,file=fname,status='unknown',iostat=ioerr)
      if (ioerr /= 0) then
         write (*,*) "Error ",ioerr," opening file ",fname
         stop
      end if
      call prthdr(10,dT1,nT1,crate0,run_title,xgridf)
      call pcnvrg(10,eQ,dQ,sQ)
      if (crate .le. crate0) then
            write (10,'("#",/,"# Global convergence achieved at T1 = ", &
                 & 1p,e8.2," : ",1p,e8.2, &
                 & " = crate < crate0 = ",e8.2)') T1,crate,crate0
      else 
          write (10,'("#",/,"# Global convergence NOT achieved by T1 = ", &
               & 1p,e8.2," : ",1p,e8.2, &
               & " = crate > crate0 = ",e8.2)') T1,crate,crate0
      end if
!dbg      write (10,'("#",1(a4,1x),1(a12,1x),28(a11,i1,1x))') & 
      write (10,'("#",1(a4,1x),1(a12,1x),193(a11,i1,1x))') & 
           '#X','X',('U',k,'UX',k,'UXX',k,k=1,nM),'Trel'
      do i = 1,NSOL
!dbg         write (10,'(1(1x,i4),1p,1(1x,e12.5),(28(1x,e12.5)))') &
         write (10,'(1(1x,i4),1p,1(1x,e12.5),(193(1x,e12.5)))') &
              i,XSOL(i),((USOL(k,i,id),id=1,NDERV+1),k=1,nM), &
         trel(i,iMtyp)
      enddo 
      write (10,'("# T1 = ",1p,e8.2)') T1
      write (10,'("# Total elapsed time = ",1p,e9.3," s")') &
           usert()-uti
      close (10,iostat=ioerr)
      if (ioerr /= 0) then
         write (*,*) "Error ",ioerr," closing file FP_last.out"
         stop
      end if
      write (*,'("# Last solution written on ",a50)') fname
!
!     Calculating the final Q on the collocation points 
!
      fname = TRIM(path)//'/Q_last.out'
      open (10,file=fname,status='unknown',iostat=ioerr)
      if (ioerr /= 0) then
         write (*,*) "Error ",ioerr," opening file ",fname
         stop
      end if
      call prthdr(10,dT1,nT1,crate0,run_title,xgridf)
      call pcnvrg(10,eQ,dQ,sQ)
      if (crate .le. crate0) then
            write (10,'("#",/,"# Global convergence achieved at T1 = ", &
                 & 1p,e8.2," : ",1p,e8.2, &
                 & " = crate < crate0 = ",e8.2)') T1,crate,crate0
      else 
          write (10,'("#",/,"# Global convergence NOT achieved by T1 = ", &
               & 1p,e8.2," : ",1p,e8.2, &
               & " = crate > crate0 = ",e8.2)') T1,crate,crate0
      end if
!dbg      write (10,'("#",1(a4,1x),1(a12,1x),36(a11,i1,1x))') &
      write (10,'("#",1(a4,1x),1(a12,1x),256(a11,i1,1x))') &
           '#X','X',('U',k,'UX',k,'UXX',k,'Q',k,k=1,nM)
!
      do i = 1,nX
!dbg         write (10,'(1x,i4,1p,1(1x,e12.5),36(1x,e12.5))') &
         write (10,'(1x,i4,1p,1(1x,e12.5),256(1x,e12.5))') &
              i,Xg(i), &
              (Ug(k,i),UXg(k,i),UXXg(k,i),Qg(k,i),k=1,nM)
      enddo 
      write (10,'("# T1 = ",1p,e8.2)') T1
      write (10,'("# Total elapsed time = ",1p,e9.3," s")') &
         usert()-uti
      close (10,iostat=ioerr)
      if (ioerr /= 0) then
         write (*,*) "Error ",ioerr," closing file Q_last.out"
         stop
      end if
      write (*,'("# Last Q solution written on ",a50)') fname
!
! Termination
!
      deallocate(WORKg1)
      deallocate(WORKg2)
      deallocate(dQ)
      deallocate(sQ)
      deallocate(eQ)
      deallocate(Ug0)
      deallocate(XSOL)
      deallocate(USOL)
      deallocate(SCTCH)
      call terminate()
      utf = usert()
      write (*,*)
      write (*,'("Program FPSolver (",a10,") ends. Elapsed time = ", &
           & 1p,e9.3)') &
           version,utf-uti
      stop
    contains
!///////////////////////////////////////////////////////////////////////
! Check convergence in loss-cone-less solution by constancy of Q,
! as estimated by two measures:
! (1) Delta Q = max_x[Q(x)]-min_x[Q(x)]
! (2) rms Q 
      subroutine Qcnvrg(Q,nM,nX,dQ,eQ,sQ)      
!                       i i    i  o  o  o
! dQ: max_x[Q(x)]-min_x[Q(x)]
! eQ: mean value of current
! sQ: rms value of current
!
        implicit none
        save	! all variables
        integer :: nM,nX
        real :: Q(nM,nX),dQ(nM),eQ(nM),sQ(nM)
!
        integer iX,iM
        real minQ,maxQ,aQ,vQ,Qx
!
        do iM = 1,nM
           minQ = +huge(minQ)
           maxQ = -huge(maxQ)
           aQ = 0.0
           vQ = 0.0         
           do iX = 1,nX
              Qx = Q(iM,iX)
              minQ = min(Qx,minQ)
              maxQ = max(Qx,maxQ)
              aQ = aQ + Qx
              vQ = vQ + Qx**2 
           enddo
           dQ(iM) = maxQ-minQ
           eQ(iM) = aQ/nX
           sQ(iM) = sqrt((vQ-aQ**2/nX)/(nX-1))
        enddo
    return 
  end subroutine Qcnvrg
!///////////////////////////////////////////////////////////////////////
! Print a standard header (as comment lines) on output unit un
  subroutine prthdr(un,dT1,nT1,crate0,run_title,xgridf)
    use FPglobal
    implicit none
    integer :: un
    real :: dT1,crate0
    integer :: nT1
    character :: run_title*100, xgridf*50
!
    integer :: k
    character :: fdate*24
!
    write (un,'("# FPSolverMS ",a10,": ",a24)') &
         version,fdate() 
!
!  Conditional inclusion of a LC term
!
#ifdef __LCTERM__ 
    write (un,'("# With a loss-cone term in the FP equation")')
#else
    write (un,'("# Without a loss-cone term in the FP equation")')
#endif
!
    write (un,'("# ",a)') run_title
    write (un,'("#")')
    write (un,'("# T0 = ",1p,e9.3," Tmax = ",e9.3, &
         & " nT = ",i3," crate0 = ",e9.3)') &
         T0,T0+dT1*nT1,nT1,crate0
    write (un,'("# x1 = ",1p,e9.3," xD = ",e9.3)') x1,xD
    write (un,'("# nM = ",i2)') nM
    do k = 1,nM
       write (un,'("#",1p," m",i1," = ",e9.3, &
            & " Cm",i1," = ",e9.3, &
            & " Dm",i1," = ",e9.3, &
            & " Sm",i1," = ",e9.3)') &
            k,m(k),k,Cm(k),k,Dm(k),k,Sm(k)
    enddo
    write (un,'("# EPS = ",1p,e9.3," DT = ",e9.3, &
         & " NINT = ",i4," NCPTS = ",i4)') &
         EPS,DT,NINT,NCPTS
    if (xgridf(1:1) .eq. ' ') then
       write (un,'("# Using automatic x-grid")')
    else 
       write (un,'("# Using given x-grid ",a50)') xgridf
    endif
    if (initf(1:1) == ' ') then
       write (un,'("# Using default initial conditions")')
    else 
       write (un,'("# Using given initial conditions ",a50)') initf
    end if
    write (un,'("# KORD = ",i2," NCC = ",i2," MF = ",i2)') KORD,NCC,MF 
    write (un,'("#")')
    return 
  end subroutine prthdr
!///////////////////////////////////////////////////////////////////////
! Print info on convergence to constant current
  subroutine pcnvrg(un,eQ,dQ,sQ)
    use FPglobal
    implicit none     
!
    integer :: un      
    real :: eQ(nM),dQ(nM),sQ(nM)
!
    integer :: k
!
    write (un,'("# <Q",i1,"> = ",1p,e9.2, &
         & " Delta Q",i1," = ",1p,e9.3, &
         & " rms_Q",i1," = ",e9.3, &
         & " : rms_Q",i1,"/<Q",i1,"> = ",e9.3)') &
         (k,eQ(k),k,dQ(k),k,sQ(k),k,k,sQ(k)/eQ(k),k=1,nM) 
    return 
  end subroutine pcnvrg
!///////////////////////////////////////////////////////////////////////
! Calculate the approximate value of the local(*) relaxation time
! using HA06b Eq. 7
! (*) Local in phase space
  real function trel(iX,iMtyp)
    use FPglobal
    implicit none
    integer :: iX,iMtyp
! 
    integer iM
!    
    trel = 0.0
    do iM = 1,nM
       trel = trel + USOL(iM,iX,1)*(m(iM)/m(iMtyp))**2
    enddo
    if (trel .gt. 0.0) then
       trel = 1./trel
    else 
       trel = huge(trel)
    endif
    return 
  end function trel
!///////////////////////////////////////////////////////////////////////
! Returns elapsed user time since start of process (sec)
  real function usert()
    implicit none
    real, save :: tarray(2),total
!
    call etime(tarray,total)
    usert = tarray(1)
    return
  end function usert
!
end program FPSolverMS
