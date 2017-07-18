!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     TO DO:
!     * FIX density normalization: n_h -> n_*(infinity) (n_h>n_*!!!)
!       Currently, the input nh is treatd as the density of the typical
!       at infinity (where phi=0), and therefore 
!       n_*,ub(rh)= "nh"*[2/sqrt(pi)+exp(1)*erfc(1)]
!       In practice, what is "known" is the TOTAL mass density at rh,
!       by kinematical measurements or the M/sig relation
!     * Read input parameters (m Cm/Sm, Dm) from FP file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!     Calculate the space density of the cusp based on the             !
!     dimensionless DF g_M(x), where x=M*/beta*(E/M)                   !
!                                                                      !
!     NOTE: at present, the population parameters (m, Cm) are not read !
!     from the FPSolver input file, but are rather typed by the user   !
!     as input to density.f. This introduces the danger of inconsis-   !
!     tencies.                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program density
      implicit none
      integer mx,nx,mm,nm
      real nh,rhoh,Nench,Mench
      integer nmode
      parameter (mx = 200,mm = 10)
      real x(mx),gm(mx,mm),gxm(mx,mm)
      real m(mm),mbh,ms,rh,Cm(mm),Dm(mm),Sm(mm)
      integer xcol,gmcol(mm),gxmcol(mm)
      character file*50,line*250
      integer mval,nval
      parameter (mval = 100)
      real val(mval)
!
      namelist /DNSINP/ 
     >   nr,nm,m,Cm,Dm,Sm,mbh,rh,ms,nmode,nh,rhoh,xcol,
     >   Nench,Mench,gmcol,gxmcol,file
!
      integer mr,nr
      parameter (mr = 100)
      real r(mr),d(mr,mm),db(mr,mm),dub(mr,mm),c(mr,mm),n(mr),rho(mr)
      real n0,ntyph,Nenc,Menc
      real dlogr,trlx,trlxh,sig2,sig2h
      integer i,im,imtyp,ioerr
      real d1,db1,dub1
      real pi
!
      real G,pc,M0,yr
      parameter (G = 6.6720e-8)
      parameter (pc = 3.08568e18)
      parameter (M0 = 1.9891e33)
      parameter (yr = 365*24*3600.)
!     
      pi = atan(1.0)*4
      do im = 1,mm
         m(im) = 0.0
         Sm(im) = -1.0
         Cm(im) = -1.0
         Dm(im) = -1.0
      enddo 
!
      nh = -1.0                 ! Total/typical number density at rh
      rhoh = -1.0               ! Total mass density at rh
      Nench = -1.0              ! Enclosed total/typical number in rh
      Mench = -1.0              ! Total enclosed stellar mass in rh
!
      read (5,nml=DNSINP)
      write (6,nml=DNSINP)
      if (nr .gt. mr) then
         print *,'ERROR: increase mr to ',nr
         stop
      endif
!
!     Either Sm and Dm are given, and Cm is calculated
!     or Cm and Dm are given, and Sm is calculated
!
      do im = 1,nm
         if (Sm(im) .lt. 0.0 .and. Cm(im) .lt. 0.0) then
            print '("Illegal data for mass ",i1)',im
            stop
         else if (Sm(im) .lt. 0.0) then
            Sm(im) = Cm(im)/Dm(im)**1.5
            print '("Setting Sm(",i1,") = ",1p,e8.2)',im,Sm(im)
         else
            Cm(im) = Sm(im)*Dm(im)**1.5
            print '("Setting Cm(",i1,") = ",1p,e8.2)',im,Cm(im)
         endif
      enddo
!     Locate the typical star in the list (Sm~1)
      imtyp = 1
      d1 = 1e30
      do im = 1,nm
         if (abs(Sm(im)-1.0) .lt. d1) then
            imtyp = im
            d1 = abs(Sm(im)-1.0)
         endif 
      enddo 
!      
      nval = xcol
      do im = 1,nm
         nval = max(nval,gmcol(im))
         nval = max(nval,gxmcol(im))
      enddo 
      if (nval .gt. mval) then
         print *,'ERROR: increase mval to ',nval
         stop
      endif 
!
      open (unit=10,file=file,status='old',iostat=ioerr,err=901)
      nx = 0
 1    continue 
         read (10,'(a250)',end=2) line
         if (line(1:1) .eq. '#') goto 1
         read (line,*,end=2) (val(i),i=1,nval)
         nx = nx+1
         if (nx .gt. mx) then
            print *,'ERROR: increase mx to at least ',nx
            stop
         endif 
         x(nx) = val(xcol)
         do im = 1,nm            
            gm(nx,im) = val(gmcol(im))
            gxm(nx,im) = val(gxmcol(im))
         enddo 
         print '(i3,1p,12(1x,e9.3))',
     >      nx,x(nx),(gm(nx,im),gxm(nx,im),im=1,nm)
         goto 1
 2    continue 
      print *,'nx = ',nx
      close (unit=10,iostat=ioerr,err=902)
!
!     It is assumed that the x values are ordered, and x(2) is the
!     first non-zero x-value
!
!     Intializing the logarithmic distance grid
!
      dlogr = exp(log(x(nx)/x(2))/nr)
      r(nr) = rh
      do i = nr-1,1,-1
         r(i) = r(i+1)/dlogr
      enddo 
!
!     Calculate the density and cumulative number
!
!     1st pass: Normalize
!
!    Note the distinction between n_*(rh) and n_*(infinity)!
!    The prefactor of the unbound DF corresponds to the TOTAL density
!    at infinity, and NOT to n_h(rh) as was erroneously assumed in
!    earlier versions of the code (and in HA06b???)
!    n_*,ub(rh)= n_*(infinity)*[2/sqrt(pi)+exp(1)*erfc(1)]
!
      do i = 1,nr
         n(i) = 0.0
         rho(i) = 0.0
         trlx = 0.0
         do im = 1,nm
            call clcnr(rh/r(i),Cm(im),Dm(im),m(im),
     >         x,gm(1,im),gxm(1,im),nx,d1,db1,dub1)
            d(i,im) = d1
            db(i,im) = db1
            dub(i,im) = dub1
            n(i) = n(i)+d(i,im)
            rho(i) = rho(i)+d(i,im)*m(im)
         enddo         
      enddo 
      do im = 1,nm
         c(1,im) = 0.0
      enddo 
      do i = 2,nr
         Nenc = 0.0
         Menc = 0.0
         do im = 1,nm
            c(i,im) = c(i-1,im) + 
     >         4*pi*((r(i)+r(i-1))/2)**2*(r(i)-r(i-1))*d(i,im)
            Nenc = Nenc+c(i,im)
            Menc = Menc+c(i,im)*m(im)
         enddo 
      enddo 
!
!     6 normalization modes:
!     nmode = 1: By given total stellar number density at rh
!     nmode = 2: By given total enclosed stellar number at rh
!     nmode = 3: By given total stellar mass density at rh
!     nmode = 4: By given total enclosed stellar mass at rh
!     nmode = 5: By given stellar number density of typical star at rh
!     nmode = 6: By given total enclosed stellar number of typical star at rh
!
      if     (nmode .eq. 1) then
         n0 = nh/n(nr)
      elseif (nmode .eq. 2) then
         n0 = Nench/Nenc
      elseif (nmode .eq. 3) then
         n0 = rhoh/rho(nr)
      elseif (nmode .eq. 4) then
         n0 = Mench/Menc
      elseif (nmode .eq. 5) then
         n0 = nh/d(nr,imtyp)
      elseif (nmode .eq. 6) then
         n0 = Nenc/c(nr,imtyp)
      else
         print *,'Illegal mormalization mode ',nmode
         stop
      endif 
!
!     The number density of the typical star at rh
!
      ntyph = n0*d(nr,imtyp)
!
!     Calculating the typical trel at rh
!
      sig2h = G*M0*mbh/(rh*pc)
!     Calculating the log to avoid numeric under/over-flow
!bug???!!!      trlxh = log(2./32/pi**2/log(mbh/ms))+1.5*log(2*pi*sig2h)-
      trlxh = log(3./32/pi**2/log(mbh/ms))+1.5*log(2*pi*sig2h)-
     >   2*log(G*M0*ms)-log(ntyph)+3*log(pc)
      trlxh = exp(trlxh)
      print '("Typical star density at rh = ",1p,e9.3," pc^-3")',ntyph
      print '("Relaxation time at rh = ",1p,e9.3," yr")',trlxh/yr      
!      
!
!     2nd pass: print out normalizaed density
!
      file = 'density.out'
      print *,'Writing output file ',file
      open (unit=10,file=file,status='unknown',iostat=ioerr,err=901)
      write (10,'(1p,"# rh = ",g9.3," pc")') rh
      write (10,'("# nmode = ",i1,1p,": nh = ",g9.3," pc^-3, rhoh = ",
     >   g9.3," Mo pc^-3, Nench = ",g9.3,", Mench = ",g9.3," Mo")') 
     >   nmode,nh,rhoh,nench,mench
      write (10,'(1p,"# m, Cm, Dm =",5(3(1x,g8.2),:,","))') 
     >   (m(im),Cm(im),Dm(im),im=1,nm)
      write (10,'("#")')
      write (10,'("#",7(a9,1x))') 
     >   'r [pc]',('n [pc^-3]',im=1,nm),'n[Mo/pc3]','trel [yr]'
      file = 'density_ub.out'
      print *,'Writing output file ',file
      open (unit=11,file=file,status='unknown',iostat=ioerr,err=901)
      write (11,'(1p,"# rh = ",g9.3," pc")') rh
      write (11,'("# nmode = ",i1,1p,": nh = ",g9.3," pc^-3, rhoh = ",
     >   g9.3," Mo pc^-3, Nench = ",g9.3,", Mench = ",g9.3," Mo")') 
     >   nmode,nh,rhoh,nench,mench
      write (11,'(1p,"# m, Cm, Dm =",5(3(1x,g8.2),:,","))') 
     >   (m(im),Cm(im),Dm(im),im=1,nm)
      write (11,'("#")')
      write (11,'("#",5(a9,1x))') 
     >   'r [pc]',('n [pc^-3]',im=1,nm)
      do i = 1,nr
         trlx = 0.0
         do im = 1,nm
            trlx = trlx+(m(im)/m(imtyp))**2*d(i,im)
            d(i,im) = n0*d(i,im)
            db(i,im) = n0*db(i,im)
            dub(i,im) = n0*dub(i,im)
         enddo         
         n(i) = n0*n(i)
         rho(i) = n0*rho(i)
!     Local relaxation time [s]
!     Scaling the relaxation time at rh to the contribution of
!     the local mass function and v.d.
         trlx = trlxh/trlx*(rh/r(i))**1.5
         write (10,'(1p,7(1x,e9.3))') 
     >      r(i),(d(i,im),im=1,nm),rho(i),trlx/yr
         write (11,'(1p,20(1x,e9.3))') 
     >      r(i),(dub(i,im),im=1,nm),(dub(i,im)/dub(i,imtyp),im=1,nm)
      enddo 
      close (unit=11,iostat=ioerr,err=902)
      close (unit=10,iostat=ioerr,err=902)
!
!     Calculating the enclosed number as function of radius
!
      file = 'number.out'
      print *,'Writing output file ',file
      open (unit=10,file=file,status='unknown',iostat=ioerr,err=901)
      write (10,'(1p,"# rh = ",g9.3," pc")') rh
      write (10,'("# nmode = ",i1,1p,": nh = ",g9.3," pc^-3, rhoh = ",
     >   g9.3," Mo pc^-3, Nench = ",g9.3,", Mench = ",g9.3," Mo")') 
     >   nmode,nh,rhoh,nench,mench
      write (10,'(1p,"# m, Cm, Dm =",5(3(1x,g8.2),:,","))') 
     >   (m(im),Cm(im),Dm(im),im=1,nm)
      write (10,'("#")')
      write (10,'("#",6(a9,1x))') 'r [pc]',('N(<r)',im=1,nm),'M*<r [Mo]'
!
      do i = 1,nr
         Menc = 0
         do im = 1,nm
            c(i,im) = n0*c(i,im)
            Menc = Menc+c(i,im)*m(im)
         enddo 
         write (10,'(1p,6(1x,e9.3))') r(i),(c(i,im),im=1,nm),Menc
      enddo 
      close (unit=10,iostat=ioerr,err=902)
!
      stop 
!
 901  continue 
      print *,'Error ',ioerr,' opening ',file
      stop
 902  continue 
      print *,'Error ',ioerr,' closing ',file
      stop
      end      
!///////////////////////////////////////////////////////////////////////
!     Calculate the density at r (up to the normalization density n*)
!     NOTE: Here xpot = rh/r, which assumes equipartition 
!     (rh of typical star). Generally, xpot = (rh*beta_*)/(r*beta_m)
!
!     dtot: total stellar number density (up to n*)
!     db  : bound stars number density (up to n*)
!     dub : unbound stars number density (up to n*)
      subroutine clcnr(xpot,Cm,Dm,m,x,gm,gxm,nx,dtot,db,dub)
!                      i  i  i i  i i  i   i i  o    o  o
      implicit none
      real xpot,Cm,Dm,m,dtot,db,dub
      integer nx
      real x(nx),gm(nx),gxm(nx)
!
      real dxint
!
      real f1,dfdx1,f2,dfdx2
      integer ix
      real z,pi
!     
      pi = atan(1.0)*4
!
!     Integrating numerically (by cubic spline) over the bound stars
!     db = (2/sqrt(pi))*I_0^xpot[gm(x)*sqrt(xpot-x)]dx
!
      db = 0.0
      do ix = 2,nx
         if (xpot .le. x(ix) .or. xpot .le. x(ix-1)) exit
         f1 = gm(ix-1)*sqrt(xpot-x(ix-1))
         dfdx1 = gxm(ix-1)*
     >      sqrt(xpot-x(ix-1))-gm(ix-1)/2/sqrt(xpot-x(ix-1))
         f2 = gm(ix)*sqrt(xpot-x(ix))
         dfdx2 =  gxm(ix)*
     >      sqrt(xpot-x(ix))-gm(ix)/2/sqrt(xpot-x(ix))
         db = db + dxint(x(ix)-x(ix-1),f1,dfdx1,f2,dfdx2)
      enddo 
      db = (2./sqrt(pi))*db
!
!     Analytic formula for the unbound stars (equipartition assumed!)
!
!     z = exp(m*xpot)*erfc(sqrt(m*xpot))
!     calculated via the log for numeric underflow protection
      z = exp(m*xpot+log(erfc(sqrt(m*xpot))))      
!bug!!      sub = 2*sqrt(pi)*Cm*(sqrt(xpot)/m+sqrt(pi)/2*z/m**1.5)
      dub = Cm/Dm**1.5*(2./sqrt(pi)*sqrt(m*xpot)+z)
!
      dtot = db+dub
      return 
      end
!///////////////////////////////////////////////////////////////////////
!     Integrate by cubic hermite spline the function f, with derivative
!     dfdx over the interval dx, where x2=x1+dx, 
!     f1=f(x1), dfdx1=f'(x1), f2=f(x2), dfdx2=f'(x2) 
      real function dxint(dx,f1,dfdx1,f2,dfdx2)
      implicit none
      real dx,f1,dfdx1,f2,dfdx2
!
      dxint = dx*((f1+f2)/2 + dx*(dfdx1-dfdx2)/12)
      return 
      end
