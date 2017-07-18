!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module FPcurrent
contains 
!///////////////////////////////////////////////////////////////////////
! Initialize variables needed by Qgrid to make computation efficient 
! Note: The nX collocation points X are in WORK(IW3)
  subroutine Qginit(X)
    use FPglobal
    use FPvars
    implicit none
    save	! all variables      
    real :: X(nX)
!
    integer :: im,iX
!
    do iX = 1,nX
       x32(iX) = sqrt(X(iX))**3
       do im = 1,nM
          Xm(im,iX) = m(im)/x32(iX)
          Ym(im,iX) = m(im)**2*Cm(im)/Dm(im)/x32(iX)
       enddo
    enddo
    return 
  end subroutine Qginit
!///////////////////////////////////////////////////////////////////////
! An efficient calculation procedure of the currents on the 
! collocation grid
  subroutine Qgrid(Q,X,g,gX,gXX,z,dzdy,I1m,I2m,I3m,RRm,SSm)
!                  o i i i  i   w w    w   w   w   w   w        
    use FPglobal
    use FPvars
    implicit none
    save	! all variables
!
    real :: X(nX),g(nM,nX),gX(nM,nX),gXX(nM,nX),Q(nM,nX)
!
! Work areas for the calculations
! Note the the index order I?m(nX,nM) is reversed wrt 
! g,gX,gXX,Q,Xm,Ym (nM,nX)
!
    real :: z(nX),dzdy(nX)
    real :: I1m(nX,nM),I2m(nX,nM),I3m(nX,nM)
    real :: RRm(nM,nX),SSm(nM,nX)
!
    integer :: im,iQM,iX
    real :: sumRRm,sumSSm,m1,m2
!
    integer, parameter :: direct = 1, revers = -1
!
! DEBUG: is this the correct place to initialize?    
    if (QCinit) then
       call  Qginit(X)
       QCinit = .false. 
    end if
!
! Calculate the auxilary integrals I1m, I2m, I3m
! Note the the index order I?m(nX,nM) is reversed wrt 
! g,gX,gXX,Q,Xm,Ym (nM,nX)
!
    do im = 1,nM
       m1 = m(im)
       m2 = m1**2
!
! Calculate the cumulative integral of 
! I1m(x) = (m^2/x^(3/2))*I[g(y)]dy(0..x)
!
       do iX = 1,nX
          z(iX) = g(im,iX)
          dzdy(iX) = gX(im,iX)
       end do
       call cubeh(nX,X,z,dzdy,I1m(1,im),direct)
! Dealing separately with X(1)=0 (kluge?!)
       I1m(1,im) = m2*I1m(1,im)
       do iX = 2,nX
          I1m(iX,im) = (m2/x32(iX))*I1m(iX,im)
       end do
!
! Calculate the cumulative integrals of 
! I2m(x) = m*I[{dg(y)/dy}/y^(3/2)]dy(x..xD)
!
! Dealing separately with y(1)=0 (kluge?!)
       z(1) = 0.0
       dzdy(1) = 0.0
       do iX = 2,nX
          z(iX) = gX(im,iX)/x32(iX)
          dzdy(iX) = (gXX(im,iX)-gX(im,iX)*1.5/X(iX))/x32(iX)
       end do
       call cubeh(nX,X,z,dzdy,I2m(1,im),revers)
       do iX = 1,nX
          I2m(iX,im) = m1*I2m(iX,im)
       end do
!    
! Calculate the cumulative integrals of 
! I3m(x) = m^2*I[g(y)/y^(3/2)]dy(x..xD)     
!     
! Dealing separately with y(1)=0 (kluge?!)
       z(1) = 0.0
       dzdy(1) = 0.0
       do iX = 2,nX
          z(iX) = g(im,iX)/x32(iX)
          dzdy(iX) = (gX(im,iX)-g(im,iX)*1.5/X(iX))/x32(iX)
       end do
       call cubeh(nX,X,z,dzdy,I3m(1,im),revers)
       do iX = 1,nX
          I3m(iX,im) = m2*I3m(iX,im)
       end do
!         
    end do                     ! End of loop on im
!
! Calculating the auxiliary sum terms Rm and Sm 
! Note the the index order I?m(nX,nM) is reversed wrt 
! g,gX,gXX,Q,Xm,Ym (nM,nX)
!
    do im = 1,nM
       do iX = 1,nX
          RRm(im,iX) = Xm(im,iX)*g(im,iX)+I2m(iX,im)
          SSm(im,iX) = Ym(im,iX)+I1m(iX,im)+I3m(iX,im)
       end do
    end do
!
! Calculating the currents
!
    do iQM = 1,nM
       do iX = 1,nX
          sumRRm = 0.0
          sumSSm = 0.0
          do im = 1,nM
             sumRRm = sumRRm + RRm(im,iX)
             sumSSm = sumSSm + SSm(im,iX)
          end do
          Q(iQM,iX) = m(iQM)*g(iQM,iX)*sumRRm-gX(iQM,iX)*sumSSm
       end do
    end do
! A kluge: set the current on the unbound boundary equal to that 
! of the first positive X
    do iQM = 1,nM
       Q(iQM,1) =  Q(iQM,2)
    end do
    return 
  end subroutine Qgrid
!///////////////////////////////////////////////////////////////////////
! An efficient calculation procedure of the currents on the 
! collocation grid
  subroutine Qgrid1(xi,Q,X,g,gX,gXX,z,dzdy,I1m,I2m,I3m,RRm,SSm)
!                   i  o i i i  i   w w    w   w   w   w   w
    use FPglobal
    use FPvars
    implicit none

    save	! all variables
!
    integer :: xi
    real :: X(nX),g(nM,nX),gX(nM,nX),gXX(nM,nX),Q(nM,nX)
!
! Work areas for the calculations
! Note the the index order I?m(nX,nM) is reversed wrt 
! g,gX,gXX,Q,Xm,Ym (nM,nX)
!
    real :: z(nX),dzdy(nX)
    real :: I1m(nX,nM),I2m(nX,nM),I3m(nX,nM)
    real :: RRm(nM,nX),SSm(nM,nX)
!
    integer :: im,iQM,iX
    real :: sumRRm,sumSSm,m1,m2
    integer :: xp1,xm1,x2
!
    integer, parameter :: direct = 1, revers = -1
!
    logical init
    data init/.true./
!
! DEBUG: is this the correct place to initialize?    
    if (init) then
       call  Qginit(X)
       init = .false. 
    end if
!
    xp1 = min(xi+1,nX)
    xm1 = max(xi-1,1)
    x2 = max(xm1,2)
!
! Calculate the auxilary integrals I1m, I2m, I3m
! Note the the index order I?m(nX,nM) is reversed wrt 
! g,gX,gXX,Q,Xm,Ym (nM,nX)
!
    do im = 1,nM
       m1 = m(im)
       m2 = m1**2
!
! Calculate the cumulative integral of 
! I1m(x) = (m^2/x^(3/2))*I[g(y)]dy(0..x)
!
       do iX = 1,xp1
          z(iX) = g(im,iX)
          dzdy(iX) = gX(im,iX)
       end do
       call cubehp(nX,1,xp1,X,z,dzdy,I1m(1,im),direct)
! Dealing separately with X(1)=0 (kluge?!)
       I1m(1,im) = m2*I1m(1,im)
       do iX = 2,xp1
          I1m(iX,im) = (m2/x32(iX))*I1m(iX,im)
       end do
!
! Calculate the cumulative integrals of 
! I2m(x) = m*I[{dg(y)/dy}/y^(3/2)]dy(x..xD)
!
! Dealing separately with y(1)=0 (kluge?!)
       z(1) = 0.0
       dzdy(1) = 0.0
       do iX = x2,nX
          z(iX) = gX(im,iX)/x32(iX)
          dzdy(iX) = (gXX(im,iX)-gX(im,iX)*1.5/X(iX))/x32(iX)
       end do
       call cubehp(nX,xm1,nX,X,z,dzdy,I2m(1,im),revers)
       do iX = xm1,nX         ! dbg xm1 or x2?
          I2m(iX,im) = m1*I2m(iX,im)
       end do
!    
! Calculate the cumulative integrals of 
! I3m(x) = m^2*I[g(y)/y^(3/2)]dy(x..xD)     
!     
! Dealing separately with y(1)=0 (kluge?!)
       z(1) = 0.0
       dzdy(1) = 0.0
       do iX = x2,nX
          z(iX) = g(im,iX)/x32(iX) ! redundant: same as for I2m!
          dzdy(iX) = (gX(im,iX)-g(im,iX)*1.5/X(iX))/x32(iX)
       end do
       call cubehp(nX,xm1,nX,X,z,dzdy,I3m(1,im),revers)
       do iX = xm1,nX         ! dbg xm1 or x2?
          I3m(iX,im) = m2*I3m(iX,im)
       end do
!         
    end do                     ! End of loop on im
!
! Calculating the auxiliary sum terms Rm and Sm 
! Note the the index order I?m(nX,nM) is reversed wrt 
! g,gX,gXX,Q,Xm,Ym (nM,nX)
!
    do im = 1,nM
! 2nd order estimate of dQ/dx requires evaluation of Q at xi+/-1
       do iX = x2,xp1,2
          RRm(im,iX) = Xm(im,iX)*g(im,iX)+I2m(iX,im)
          SSm(im,iX) = Ym(im,iX)+I1m(iX,im)+I3m(iX,im)
       end do
    end do
!
! Calculating the currents
!
    do iQM = 1,nM
! 2nd order estimate of dQ/dx requires evaluation of Q at xi+/-1
       do iX = x2,xp1,2
          sumRRm = 0.0
          sumSSm = 0.0
          do im = 1,nM
             sumRRm = sumRRm + RRm(im,iX) 
             sumSSm = sumSSm + SSm(im,iX)
          end do
          Q(iQM,iX) = m(iQM)*g(iQM,iX)*sumRRm-gX(iQM,iX)*sumSSm
       end do
    end do
! A kluge: set the current on the unbound boundary equal to that 
! of the first positive X
    do iQM = 1,nM
       Q(iQM,1) =  Q(iQM,2)
    end do
    return 
  end subroutine Qgrid1
!///////////////////////////////////////////////////////////////////////
! Construct the cumulative function c(x_i) for a function f(x_i)
! given on an equally-spaced grid {x_i}, with its 1st derivative
! f'(x_i), using the cubic hermite spline.
! The grid x is assumed to be ordered x_i<x_(i+1)
!
! see http://en.wikipedia.org/wiki/Catmull-Rom_spline
! and mathematica notebook ./cubeh.nb
!
  subroutine cubeh(n,x,f,dfdx,c,direct)
!                  i i i i    o i
    implicit none
    integer, intent(in) :: n,direct
    real, intent(in) :: x(n),f(n),dfdx(n)
    real, intent(out) :: c(n)
!
    integer :: i
!
    if (direct .eq. 1) then
!     Calculate the direct cumulative function c(<x)
       c(1) = 0.0
       do i = 2,n,+1
          c(i) = c(i-1)+ &
               dxint(x(i)-x(i-1),f(i-1),dfdx(i-1),f(i),dfdx(i))
       enddo
    else 
!     Calculate the inverse cumulative function c(>x)
       c(n) = 0.0
       do i = n-1,1,-1
          c(i) = c(i+1)+ &
               dxint(x(i+1)-x(i),f(i),dfdx(i),f(i+1),dfdx(i+1))
       enddo
    endif
    return 
  end subroutine cubeh
!///////////////////////////////////////////////////////////////////////
! Construct the cumulative function c(x_i) for a function f(x_i)
! given on an equally-spaced grid {x_i}, with its 1st derivative
! f'(x_i), using the cubic hermite spline.
! The grid x is assumed to be ordered x_i<x_(i+1)
!
! see http://en.wikipedia.org/wiki/Catmull-Rom_spline
! and mathematica notebook ./cubeh.nb
!
! This version of cubeh can calculate partial cumulative functions
! from x(x1) to x(xn) (direct) or from x(xn) to x(x1) (reverse) 
! to save computing time
!
  subroutine cubehp(n,x1,xn,x,f,dfdx,c,direct)
!                   i i  i  i i i    o i
    implicit none
    integer, intent(in) :: n,x1,xn,direct
    real, intent(in) :: x(n),f(n),dfdx(n)
    real, intent(out) :: c(n)
!
    integer :: i
!
    if (direct .eq. 1) then
! Calculate the direct cumulative function c(<x)
       c(1) = 0.0
       do i = 2,xn,+1
          c(i) = c(i-1)+ &
               dxint(x(i)-x(i-1),f(i-1),dfdx(i-1),f(i),dfdx(i))
       enddo
    else 
! Calculate the inverse cumulative function c(>x)
       c(n) = 0.0
       do i = n-1,x1,-1
          c(i) = c(i+1)+ &
               dxint(x(i+1)-x(i),f(i),dfdx(i),f(i+1),dfdx(i+1))
       enddo
    endif
    return 
  end subroutine cubehp
!///////////////////////////////////////////////////////////////////////
! Integrate by cubic hermite spline the function f, with derivative
! dfdx over the interval dx, where x2=x1+dx, 
! f1=f(x1), dfdx1=f'(x1), f2=f(x2), dfdx2=f'(x2) 
  real function dxint(dx,f1,dfdx1,f2,dfdx2)
    implicit none
    real :: dx,f1,dfdx1,f2,dfdx2
!
    dxint = dx*((f1+f2)/2 + dx*(dfdx1-dfdx2)/12)
    return 
  end function dxint
!
end module FPcurrent
