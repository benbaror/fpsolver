! User supplied routines for PDECOL
! SUBROUTINE UINIT: Initial conditions
! SUBROUTINE BNDRY: Boundary conditions
! FUNCTION LCterm:  Loss-cone term 
MODULE FPuser
contains
!///////////////////////////////////////////////////////////////////////
! A new interface for calculating the RHS of the PDE,
! by getting U(_xi), UX(x_i), UXX(x_i) for the entire grid x_i,
! and calculating the RHS at all x_i.
! However, if xi > 0, it is assumed that only results for the xi'th 
! point are needed, and the calculations will be done more 
! economically.
  subroutine FGRID(xi,T)
!                  i  i
    use FPglobal
    use FPvars
    use FPgrid
    use FPcurrent
    implicit none
    save ! all variables
!
    integer,intent(in) :: xi
    real, intent(in) :: T
!
    integer :: iM
!
! Work areas for the current calculations
! NOTE! such work areas are also defined in the main. 
! Can they be used here for economy???
!
    REAL :: WORKg1(nX,3),WORKg2(nX*nM,5),I2m(nX,nM)
!
    integer :: i,i1,i2
    logical :: fullg
!
!     Calculating the currents
!
    fullg = xi .le. 0
    if (fullg) then
       call Qgrid(Qg,Xg,Ug,UXg,UXXg, &
            WORKg1(1,1),WORKg1(1,2), &
            WORKg2(1,1),WORKg2(1,2),WORKg2(1,3),WORKg2(1,4),WORKg2(1,5))
    else
       call Qgrid1(xi,Qg,Xg,Ug,UXg,UXXg, &
            WORKg1(1,1),WORKg1(1,2), &
            WORKg2(1,1),WORKg2(1,2),WORKg2(1,3),WORKg2(1,4),WORKg2(1,5))
    END IF
    I2m = RESHAPE(WORKg2(:,2),(/nX,nM/))

!
! DEBUG NOTE: try also direct formula for dU/dt
!
! Calculating dU(x_i)/dt = -x_i^(5/2)*dQ(x_i)/dx for all x_i
! by 2nd-order finite differencing
!
      
    if (fullg) then
       i1 = 2
       i2 = nX-1
    else if (xi .gt. 0 .and. xi .le. nX) then
       i1 = max(xi,2)
       i2 = min(xi,nX-1)
    else 
       print *,'FGRID: xi = ',xi,' > nX'
       stop
    end if
    do iM = 1,nM
       UTg(iM,1) = 0.0
       do i = i1,i2
!
! Conditional inclusion of a LC term
!
#ifdef __LCTERM__ 
          UTg(iM,i) = &
               -sqrt(Xg(i))**5* &
               (Qg(iM,i+1)-Qg(iM,i-1))/(Xg(i+1)-Xg(i-1)) &
               -LCterm(iM,i,I2m)
#else
          UTg(iM,i) = &
               -sqrt(Xg(i))**5* &
               (Qg(iM,i+1)-Qg(iM,i-1))/(Xg(i+1)-Xg(i-1))
#endif
       end do
       UTg(iM,nX) = 0.0
    end do                     ! End of loop on iM
!
!dbg(
    if (DEBUG) then                                            !dbg!
       print '("Xg   = ",1p,5(1x,e12.5)," ... ",5(1x,e12.5))', & !dbg!
            (Xg(i),i=1,5),(Xg(i),i=nX-4,nX)              !dbg!
       print '("Ug   = ",1p,5(1x,e12.5)," ... ",5(1x,e12.5))', & !dbg!
            (Ug(1,i),i=1,5),(Ug(1,i),i=nX-4,nX)          !dbg!
       print '("UXg  = ",1p,5(1x,e12.5)," ... ",5(1x,e12.5))', & !dbg!
            (UXg(1,i),i=1,5),(UXg(1,i),i=nX-4,nX)        !dbg!
       print '("UXXg = ",1p,5(1x,e12.5)," ... ",5(1x,e12.5))', & !dbg!
            (UXXg(1,i),i=1,5),(UXXg(1,i),i=nX-4,nX)      !dbg!
       print '("Qg   = ",1p,5(1x,e12.5)," ... ",5(1x,e12.5))', & !dbg!
            (Qg(1,i),i=1,5),(Qg(1,i),i=nX-4,nX)          !dbg!
       print '("UTg  = ",1p,5(1x,e12.5)," ... ",5(1x,e12.5))', & !dbg!
            (UTg(1,i),i=1,5),(UTg(1,i),i=nX-4,nX)        !dbg!
    end if                                                      !dbg!
!dbg)
    return
  end subroutine FGRID
!-----------------------------------------------------------------------
! A subroutine masquarading as the original F subroutine 
! (for calculating the LHS of the evolution equation at one
! collocation grid point iXg)
!     
! It is assumed that the grid was last initialized through the
! main subroutine.
  SUBROUTINE F(iXg0,FVAL)
!              i      o
    use FPglobal
    use FPvars
    use FPgrid
    implicit none
    save	! all variables
!     
    integer iXg0
    real FVAL(nM)
!
    integer :: iM
!     
    do iM = 1,nM
       FVAL(iM) = UTg(iM,iXg0)
    end do
    return
  end SUBROUTINE F
!///////////////////////////////////////////////////////////////////////
! Setting the boundary conditions on the left and right
!
! DEPENDING ON THE TYPE OF PDE(S), 0, 1, OR 2 BOUNDARY CONDITIONS
! ARE REQUIRED FOR EACH PDE IN THE SYSTEM. THESE ARE IMPOSED AT XLEFT
! AND/OR XRIGHT AND EACH MUST BE OF THE FORM....
!
!     B(U,UX)  =  Z(T)
!
! WHERE  B  AND  Z  ARE ARBITRARY VECTOR VALUED FUNCTIONS WITH
! nM COMPONENTS AND  U, UX, AND T  ARE AS ABOVE.  THESE BOUNDARY
! CONDITIONS MUST BE CONSISTENT WITH THE INITIAL CONDITIONS.
!

  SUBROUTINE BNDRY(T, X, U, UX, DBDU, DBDUX, DZDT)
    use FPglobal
    use FPvars
    implicit none
    save	! all variables
    real :: T,X
    real :: U(nM),UX(nM),DZDT(nM),DBDU(nM,nM),DBDUX(nM,nM)      

! NOTE: THE INCOMING VALUE OF X WILL BE EITHER XLEFT OR XRIGHT.
! IF NO BOUNDARY CONDITION IS DESIRED (NULL BOUNDARY) FOR SAY 
! THE K-TH PDE AT ONE OR BOTH OF THE ENDPOINTS XLEFT OR XRIGHT, 
! THEN DBDU(K,K) AND DBDUX(K,K) SHOULD BOTH BE SET TO ZERO.
!
    integer :: iM
!
    do iM = 1,npde
       if     (x .eq. xU) then
! U(xU) = Cm*exp(xU) = Cm
          DBDU(iM,iM) = Cm(iM)
          DBDUX(iM,iM) = 0.0
          DZDT(iM) = 0.0
       else if (x .eq. xD) then
! U(xR) = 0
          DBDU(iM,iM) = 1.0
          DBDUX(iM,iM) = 0.0
          DZDT(iM) = 0.0
       else 
          write (*,*) "FATAL ERROR IN BNDRY:"
          write (*,*) "X = ",X," <> XU = ",xU," <> XD = ",xD
          stop
       endif
    end do
    return 
  end SUBROUTINE BNDRY
!///////////////////////////////////////////////////////////////////////
! Initial configuration, consistent with the boundary conditions 
  SUBROUTINE UINIT(X,U)
    use FPglobal
    use FPvars
    use FPgrid
    implicit none
    save	! all variables
    real :: X,U(nM)
!
    integer :: iM
!
! Areas for the initial grid option
!
    integer :: i,ioerr
    character(len=500) :: line
    logical :: init = .true. 
    integer, parameter :: max_val = 1+9*4
    real :: val(max_val)
    integer, parameter :: max_ini = 200
    real :: x_ini(max_ini),g_ini(nM,max_ini)
    integer :: n_ini
    real w
!
! An arbitrary initial configuration which is 1 at xU and 0 at xD
!
    if (initf(1:1) == ' ') then
       do iM = 1,nM
          if     (X .eq. xU) then
             U(iM) = Cm(iM)
          else if (X == xD) then
             U(iM) = 0.0
          else 
             U(iM) = Cm(iM)*(1.0-(X-xU)/(xD-xU))
          end if
       end do
!
! Read initial conditions from a file (line starting with # are comments)
! It is assumed that the values are ordered by increasing x, and
! that the requested values of x are within the grid limit
!
    else
       if (init) then
          init = .false. 
          open (15,file=initf,status='old',iostat=ioerr)
          if (ioerr /= 0) then
             write (*,*) "Error ",ioerr," opening init file ",initf
          end if
          n_ini = 0
          do
             read(15,'(a500)',iostat=ioerr) line
             if(ioerr > 0) then
                write (*,*) "Error ",ioerr," reading init file ",initf
             else if (ioerr < 0) then
                exit
             end if
             if (line(1:1) == '#') cycle
             n_ini = n_ini+1
             read (line,*) (val(i),i=1,g_idx(nM))
             x_ini(n_ini) = val(x_idx)
             do iM = 1,nM
                g_ini(iM,n_ini) = val(g_idx(iM))
             end do
          end do
          close (15,iostat=ioerr)
          if (ioerr /= 0) then
             write (*,*) "Error ",ioerr," closing init file ",initf
          end if
          write (*,*) "Using initial solution from file ",initf
!dbg(
!          do i = 1,n_ini
!             write (*,*) i,x_ini(i),(g_ini(iM,i),iM=1,nM)
!          end do
!dbg)
       end if
!
       do iM = 1,nM
          U(iM) = -1.0 ! <0 values are illegal
          do i = 1,n_ini-1
             if (X >= x_ini(i) .and. X <= x_ini(i+1)) then
! A linear interpolation
                w = (X-x_ini(i))/(x_ini(i+1)-x_ini(i))
                U(iM) = g_ini(iM,i)*(1-w)+g_ini(iM,i+1)*w
                exit
             end if
          end do
       end do
       if (U(iM) < 0.0) then
          write (*,*) "Error in UINIT: X value ",X," out of bounds"
          stop
       end if
    end if
    return 
  end SUBROUTINE UINIT
!
! Conditional compilation of the LC term
!
#ifdef __LCTERM__ 
!///////////////////////////////////////////////////////////////////////
! Calculate the approximate non-resonant LC term for direct infall
! through the event horizon of a MBH (HA06, eqs 6 & 7)
  REAL FUNCTION LCterm(iM,iX,I2m)
    use FPglobal
    use FPvars
    use FPgrid
    implicit none 
!
    INTEGER, INTENT(in) :: iM,iX
    REAL, DIMENSION(nX,nM), INTENT(IN) :: I2m
!
    INTEGER :: jM, i,iXh(1)
    REAL :: trlx,tRR,tM, tGR_inv, Jc2, alpha, q_ha,Nm2,Nm1,Nmh
    !REAL, DIMENSION (10000) :: ecc2, tGR
    
!
    REAL, PARAMETER :: X_pinhole = 0.5


! Calculate the local relaxation time

!

!
! LC term = R_m = g_m(x)/trlx/log(J_c(x)/J_lso)
! Jc/Jlso = 2^(-5/2)*(c/sigma_*)/sqrt(x) = sqrt(r_h/r_g)/sqrt(32*x)
! where r_h = GM/(sigma_*)^2 and r_g = GM/c^2
!    
    trlx = 0.0
    DO jM = 1,nM
       trlx = trlx + Ug(jM,iX)*m(jM)**2
    ENDDO
    trlx = m(iMtyp)**2/trlx*ANR/Aeps

#ifdef pinhole
    LCterm = Ug(iM,iX)/trlx/LOG(SQRT(rh2rg/32/Xg(iX)))
    IF (Xg(iX) <= X_pinhole) THEN
      LCterm = 0.0
      !return 
   ENDIF
#else    ! BW77
   alpha = 4e6**(-2)*4e4*1.9**3*xD
   q_ha = 1.6*LOG(6.0/3.14*4e6)*alpha*Xg(iX)**(-5./2.)*Ug(iM,iX)
   LCterm =  ABS(0.5*Ug(iM,iX)*Xg(iX)**(5./2.)/alpha/LOG(4e6)/(5.56&
        &+LOG(xD/4./Xg(iX))/q_ha))
#endif pinhole

   ! Binary capture
   if (Xg(iX) < xmax .and. Xg(iX) > xmin) then
      LCterm =   - eta_L*Xg(iX)**(3./2.)/log(xmax/xmin)*sqrt(2.)/3/pi
   end if
   
#ifdef do_RR
    tM = 0.0
    Nm2 = 0.0
    Nm1 = 0.0
    Nmh = 0.0
    iXh = MINLOC(ABS(Xg - 0.5))
    DO jM = 1,nM
    !   tM = tM + Ug(jM,iX)*m(jM)
       Nm2 = Nm2 + 2/3.*(Xg(iX)**(-3/2.)*Ug(jM,iX)*m(jM)**2 + I2m(iX,jM)&
            &*m(jM))
       Nm1 = Nm1 + 2/3.*(Xg(iX)**(-3/2.)*Ug(jM,iX)*m(jM) + I2m(iX,jM))
       Nmh = Nmh + 2/3.*(Xg(iXh(1))**(-3/2.)*Ug(jM,iXh(1))*m(jM) + I2m(iXh(1),jM))
       !Ng = Ng + I2m(jM)/m(jM)
    ENDDO
!    WRITE (*,*) Xg(iX),Ng,Ug(iM,iX)
    !tM = SQRT(32.0)*SUM(Cm*m)/tM*Xg(iX)**(3.0/2.0)
    tM = Nmh/Nm1
    Jc2 = rh2rg/32./Xg(iX)
    tGR_inv = 3./8.*LOG(Jc2)/Jc2/AGR
    !trr = (1/tm+tGR_inv)*trlx*ARR/ANR*LOG(4e6)
    trr = ARR*16/3.0*LOG(4e6)*(1/tm+tGR_inv)/Nm2/Xg(ix)**(3&
         &/2.)
!    WRITE(*,*) Xg(ix),trlx,trr
    LCterm = LCterm + Ug(iM,iX)/trr
#endif
    return 
  end function LCterm
!///////////////////////////////////////////////////////////////////////
#endif
!
! Dummy definition. Function not in use.
  subroutine DERIVF(T, X, U, UX, UXX, DFDU, DFDUX, DFDUXX, nM)
    implicit none
    save	! all variables
    real :: T,X
    integer :: nM
    real :: U(nM), UX(nM), UXX(nM)
    real :: DFDU(nM,nM), DFDUX(nM,nM), DFDUXX(nM,nM)
    return 
  end subroutine DERIVF
end module FPuser
