  ! Integrate orbits in a specified potential to generate a poicare plot
  ! Version that calculates the orbit in cartesian coordinates.
  ! And uses a file read in to specify the potential and E-field.
module orbcartpoin
  real :: Bsqpsi=.85, Eropsi=.01, psi=1., Omega=0.
  real :: r0,rc0,rg0,rpmin,rpmax,rpave,wpf=0.5,B,w0=1.,aleph=1.,wpt
  real :: wpm=0.,psim=0.5   ! Minimum confined energy, for psimax
  real :: dtxb=0.04363     ! dt x Omega=2pi/(steps per cyclotron period=144)
  real :: yrfmax,zrfmax
  integer, parameter :: npmmax=100,nerpmax=20
  real, dimension(npmmax,nerpmax) :: psiplot,riplot
  real, dimension(npmmax,nerpmax,2) :: omegaplot
  integer :: nerp=1,npm=npmmax   ! n Eropsi's, psi divisions of psimax
  logical :: lp=.true.,lprint=.true.,lo=.true.,lpi=.false.
  logical :: lgread=.false.,lhread=.false.
  integer, parameter :: nbouncemax=200  ! Max bounces to store.
  integer, parameter :: nplot=25000     ! Max steps to store.
  integer, parameter :: nvec=6
  integer :: iwritetype=0,idoplots=2,irktype=0,idebug=0
  integer :: nwp=39,nstep=1000000       ! Default maximum steps to take.
  real, dimension(nplot,2) :: r,th,z,vr,vt,vz,w,tv,pt
  real, dimension(nplot,2) :: xplot,yplot,zplot,phip,Ep,phir0
  real, dimension(nplot) :: xptemp,yptemp
  real, dimension(nbouncemax) :: wp,xi,tc
  real, dimension(nvec) :: y0,y1,y2
  real ::  vangle0=0.          ! Angle of initial v to r-direction degrees.
  character*30 :: filename='helmphiguard.dat'
  data y0/10.,0.,0.,0.,0.,0./ !Initial orbit (position) defaults
contains
!***********************************************************************
   subroutine RKFUN(yp,y,t,icase,ierr)
     real yp(6),y(6),t
     integer icase,ierr
! Lorentz RK advance function. Evaluates derivative of y into yp.
     real E(3)  ! Er,Et(=0),Ez
! Evaluate E.
     call getfield(y,t,icase,E)
! Enter derivatives.
! Advance cartesian x, y, z, vr, vy
        yp(1)=y(4)        ! vr
        yp(2)=y(5)        ! vp
        yp(3)=y(6)        ! vz
        yp(4)=-E(1)-y(5)*B    ! -Er-vt*B
        yp(5)=-E(2)+y(4)*B    ! -Ep+vr*B 
        yp(6)=-E(3)
        ierr=0
    end subroutine RKFUN
!***********************************************************************
    subroutine getfield(y,t,icase,E)
!      implicit none
      real y(6),t,E(3)
      integer icase
      integer, parameter :: Lzg=100         ! Size max of read in phi
      real, dimension(-1:Lzg+1,-1:Lzg+1) :: phiguarded
      real, dimension(0:Lzg) :: zg,yg
      integer nzf,nyf
      real :: zrf,yrf,zs,rh,ys,zsa,gy,gz
      save
      if(.not.lgread)then
         call inputphi(Lzg,nzf,nyf,zg,yg,phiguarded)
         zrfmax=zg(nzf)
         zrf=(zg(nzf)-zg(0))
         yrfmax=yg(nyf)
         yrf=(yg(nyf)-yg(0))
         lgread=.true.
      endif

! Global cartesian. y(1:3) is x,y,z... So calculate z,r
      zs=(y(3)-zg(0))/zrf
      zsa=abs(zs)
      rh=sqrt(y(1)**2+y(2)**2+1.e-24)
      rs=(rh-yg(0))/yrf
! Interpolate within grid
      call grad2dint(phiguarded,nzf,nzf,nyf,min(1.,zsa),min(1.,rs),gz,gy)
! Maybe extrapolate
      if(zsa.gt.1.or.rs.gt.1)then
         earg=exp(-sqrt((max(zsa,1.)-1.)*zrf**2+((max(rs,1.)-1.)*yrf)**2))
         gz=gz*earg
         gy=gy*earg
      endif
! Project the radial electric field found into cartesian components.      
      E(1)=-gy/yrf*y(1)/rh
      E(2)=-gy/yrf*y(2)/rh
      E(3)=sign(gz/zrf,zs)         
    end subroutine getfield
!***********************************************************************
    real function phiofrz(rh,zh)
      real rh,zh
      integer, parameter :: Lzg=100         ! Size max of read in phi
      integer, parameter :: nclv=24
      real, dimension(-1:Lzg+1,-1:Lzg+1) :: phiguarded,work
      real, dimension(0:Lzg) :: zg,yg
      real, dimension(nclv) :: zclv
      real :: zrf,yrf
      external bilin2d
      save
      if(.not.lhread)then
         ! This packs as phiguarded(-1:nzf+1,-1,nyf+1)
         call inputphi(Lzg,nzf,nyf,zg,yg,phiguarded)
         zrfmax=zg(nzf)
         zrf=(zg(nzf)-zg(0))
         yrfmax=yg(nyf)
         yrf=(yg(nyf)-yg(0))
         lhread=.true.
         if(lpi)call plotinputphi(Lzg,nzf,nyf,zg,yg,nclv,zclv,phiguarded,work)
      endif
      zs=(zh-zg(0))/zrf
      zsa=abs(zs)
      ys=(rh-yg(0))/yrf
!Interpolate within grid limits
      phiofrz=bilin2d(phiguarded,nzf,nzf,nyf,min(1.,zsa),min(1.,ys))
! Possibly project beyond grid with exponential decay.
      if(zsa.gt.1.or.ys.gt.1)then
         phiofrz=phiofrz*exp(-sqrt(((max(1.,zsa)-1.)*zrf)**2 &
              +((max(1.,ys)-1.)*yrf)**2))
      endif
    end function phiofrz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine inputphi(Lz,nzf,nyf,z,y,phiguarded)
! This reads into phiguarded as if it were 
!    real phiguarded(-1:nzf+1,-1:nyf+1)
! Subsequently access phiguarded accordingly. Overriding external def.
    integer Lz,nzf,nyf
    real phiguarded(*)
    real z(0:Lz),y(0:Lz)
    open(15,file=filename,status='old',form='unformatted',err=1)
    read(15)nzf,nyf
    if(nzf.gt.Lz.or.(nzf+3)*(nyf+3).gt.(Lz+3)**2) &
        stop 'Incompatible phiguarded ranks'
    read(15)(z(i),i=0,nzf),(y(i),i=0,nyf)
    read(15)((phiguarded(i+(j-1)*(nzf+3)),i=1,nzf+3),j=1,nyf+3)
    close(15)
    write(*,'(3a,i4,f7.3,i4,f7.3)')'Read in ',filename(1:lentrim(filename)), &
         ' successfully. z, r: ',nzf,z(nzf),nyf,y(nyf)
    return
1   write(*,*)'No file ',filename(1:lentrim(filename)),   &
         ' opened. Potential read unsuccessful'
    stop
  end subroutine inputphi
!***********************************************************************
  subroutine plotinputphi(Lzg,nzf,nyf,zg,yg,nclv,zclv,phiguarded,work)
    implicit none
    integer :: Lzg,nzf,nyf
    real, dimension(-1:Lzg+1,-1:Lzg+1) :: phiguarded,work
    real, dimension(0:Lzg) :: zg,yg
    integer :: nclv,iclv,icsw,ilmax
    real, dimension(nclv) :: zclv
    ! Contour the read-in potential for verification.
    call charsize(0.018,0.018)
    call pltinaspect(0.,zg(nzf),0.,yg(nyf))
    call axis; call axis2; call axlabels('z','r')
    ilmax=nint(alog10(phiguarded(0,0)))
    do iclv=1,nclv
       zclv(iclv)=10.**(ilmax+1)*10**(-float((iclv-1)/5)) &
            *(1.-0.2*mod(iclv+4,5))
       if(iclv.eq.nclv)exit
    enddo
    icsw=1
    call contourl(phiguarded(0,0),work,Lzg+3,nzf+1,nyf+1, &
         zclv,iclv,zg,yg,icsw)
    call pltend
    call multiframe(0,0,0)  ! Cancel frame aspect.
  end subroutine plotinputphi
!***********************************************************************
!***********************************************************************
    subroutine doplots2(imax,wp0,w0,pt0)
      integer, dimension(2) :: imax
      integer, parameter :: nzplot=100
      real, dimension(nzplot) :: zlplot,psilplot,psicplot
      if(lp)then; call pfset(3); else; call pfset(-3);
      endif

      np=min(nwp,2)
      call multiframe(3,1,0)
      zmin=-8.
      zmax=8.
      call minmax(r(:,np),imax,rmin,rmax)
      rmin=rmin-.3
      rmax=rmax+.3
      call minmax(vz(:,np)**2/2.-phip(:,np),imax(np),wmin,wmax)

      nbb=nint(3.1415926/dtxb)
      call boxcarave(imax,nbb,z,xptemp)
      call boxcarave(imax,nbb,r,yptemp)
      call pltinit(zmin,zmax,rmin,rmax)
      call ticnumset(4)
      call charsize(.015,.015)
      call axis

      call winset(.true.)
      call dashset(1)
      call polyline([zmin,zmax],[rc0,rc0],2)
      call jdrwstr(wx2nx(0.9*zmin+0.1*zmax),wy2ny(rc0)-.01,'r!dc0!d',0.)
      call dashset(0)
      call color(4)
      call polyline(z,r,imax)
      call color(15)
      call polymark(z,rpmax,1,2)
      call polymark(z,rpmin,1,3)
      call winset(.false.)
      call color(2)
      call polyline(xptemp(nbb),yptemp(nbb),imax-2*nbb)
      call jdrwstr(.03,wy2ny((rmax+rmin)/2)+.05,'!p!o-!o!qr',0.1)
      call color(4);call jdrwstr(.03,wy2ny((rmax+rmin)/2),'r',-1.)
      call color(15)
      call polymark(z,r,1,1)

      nbb=nint(3.1415926/dtxb)
      call boxcarave(imax,nbb,z,xptemp)
      call boxcarave(imax,nbb,vz**2/2.-phip,yptemp)

      call pltinit(zmin,zmax,wmin,0.02)
      call charsize(.015,.015)
      call axis
      call color(4);call jdrwstr(.03,wy2ny(0.01),'W!d!A|!@!d',-1.)
      call color(2)
      call jdrwstr(.03,wy2ny((wmin+.02)*1.),'!p!o-!o!qW!d!A|!@!d',-0.1)
      call color(15)
      call winset(.true.)
      call color(4); call polyline(z,vz**2/2.-phip,imax)
      if(nwp.eq.2)then
         call color(1);
         call polyline(z(:,2),vz(:,2)**2/2.-phip(:,2),imax)
         call polymark(z(:,2),vz(:,2)**2/2.-phip(:,2),1,1)
      endif
      call color(2)
      call polyline(xptemp(nbb),yptemp(nbb),imax-2*nbb)
      call color(15)
      call polymark(z,vz**2/2.-phip,1,1)
      call dashset(4)
      call polyline([zmin,zmax],[0.,0.],2)
      call dashset(0)

      do i=1,nzplot
         zlplot(i)=zmin+(zmax-zmin)*(i-1.)/(nzplot-1)
         psicplot(i)=-phiofrz(rmin,zlplot(i))     ! At lowest r.
!         psilplot(i)=-phiofrz(r0,zlplot(i))      ! At initial r.
         psilplot(i)=-phiofrz(rc0,zlplot(i))     ! At initial gyrocenter.
      enddo
      call pltinit(zmin,zmax,min(-psi-.2*psi,psicplot(nzplot/2)), &
           -psicplot(nzplot/2)*.2)
      call charsize(.015,.015)
      call axis
      call axlabels('z','')
      call color(4); call jdrwstr(.03,wy2ny(0.),'W!d!A|!@!d',-1.)
      
      call winset(.true.)
      call color(4); call polyline(z,vz**2/2.-phip,imax)
      call color(1)
      call polyline(zlplot,psilplot,nzplot)
      call jdrwstr(.12,wy2ny(0.)-.02,'-!Af!@(r!dc0!d)',1.)
      call color(3)
      call polyline(zlplot,psicplot,nzplot)
      call jdrwstr(.12,wy2ny(0.)-.05,'-!Af!@(r!dmin!d)',1.)
!      call color(1);
!      call polyline(z(:,2),vz(:,2)**2/2.-phip(:,2),imax)
!      call polymark(z(:,2),vz(:,2)**2/2.-phip(:,2),1,1)
!      call polymark(z,vz**2/2.-phip,1,1)
      call color(15)
      call polymark(z,vz**2/2.-phip,1,1)
      call pltend

    end subroutine doplots2
!***********************************************************************
    subroutine doplots(imax,wp0,w0,pt0)
      integer, dimension(2) :: imax
      if(lp)then; call pfset(3); else; call pfset(-3);
      endif
      
      nbb=nint(3.1415926/dtxb)
      call boxcarave(imax,nbb,z,xptemp)
      call boxcarave(imax,nbb,vz**2/2.-phip,yptemp)

      call multiframe(3,1,0)
      call autoplot(z,r,imax)
      call axlabels('','r')
!      call pltend
      
      call autoplot(z,vz**2/2.-phip,imax)
      call axlabels('','v!dz!d!u2!u/2-!Af!@')
      call color(2)
      call polyline(xptemp(nbb),yptemp(nbb),imax-2*nbb);call color(15)
!      call pltend

      call autoplot(z,phir0,imax)
      call axlabels('z','!Af!@(r!dave!d)')
      call pltend

      call multiframe(4,1,0)
      call autoplot(tv,vz**2/2.-phip,imax)
      call axlabels('t','v!dz!d!u2!u/2-!Af!@')
!      call pltend

      call autoplot(tv,vz**2/2.-phir0,imax)
      call axlabels('t','v!dz!d!u2!u/2-!Af!@(r!dave!d)')
!      call pltend
      
      call autoplot(tv,(w-w(2,1)),imax)
      call axlabels('t','Dw')

      call autoplot(tv,(pt-pt0)/rc0,imax)
      call axlabels('t','Dpt')
      call pltend

      call multiframe(0,0,0)

      call charsize(.018,.018)
!      call autoplot(xplot,yplot,imax)
      call minmax(xplot,imax,xpmin,xpmax)
      call fitrange(xpmin,xpmax,6,ipow,fac10,delta,xfirst,xlast)
      call minmax(yplot,imax,ypmin,ypmax)
      call fitrange(ypmin,ypmax,6,ipow,fac10,delta,yfirst,ylast)
      call pltinaspect(xfirst,xlast,yfirst,ylast)
      call axis; call axis2
      call polyline(xplot,yplot,imax)
      call axlabels('x','y')
      call color(4)
      call polymark(xplot,yplot,1,1)
      call color(15)
      call pltend
      call multiframe(0,0,0)

      call minmax(xplot,imax,xmin,xmax)
      call minmax(yplot,imax,ymin,ymax)
      call minmax(z,imax,zmin,zmax)

1     continue
      call geteye(x2,y2,z2)
      call pltinit(0.,1.,0.,1.)
!      call scale3(xmin,xmax,ymin,ymax,-5.,5.)
      call scale3(xmin,xmax,ymin,ymax,zmin,zmax)
      call trn32(0.,0.,0.,x2,y2,z2,1)
      call axproj(igetcorner())
      call ax3labels('x','y','z')
      call poly3line(xplot,yplot,z,imax)
      isw=irotatezoom()
      call prtend('')
      if(isw.ne.0.and.isw.ne.ichar('q'))goto 1
    end subroutine doplots
end module orbcartpoin


!***********************************************************************
!****************************************************************************
   SUBROUTINE RKADVC(DX,X0,ND,YR1,icase,YR2,cond,IERR)
     use orbcartpoin
! 
! A fortran routine for ADVANCING an initial-value ODE integration
! by a fourth-order Runge-Kutta step. 
! Based upon NR but (1) a case switch, and (2) reporting conditioning.
!
! The equation to be solved is to be given in the form
!   	dY/dX = F(Y,X)
! where Y is a vector of dimension equal to the order of the system.
!    INPUTS:
!    DX	The point spacing in the independent var.
!    X0	The initial point.
!    ND	The order of the system. Must be <= 20
!    YR1	The initial vector, dimension ND (<=20)
!       icase   The case to be used in the function.
!  Outputs
!    YR2	The advanced vector, dimension ND.
!       COND    The condition of the advance defined as 
!                The magnitude difference in y-change 
!                divided by magnitude y-change, all squared.
!       ierr
!
!   RKFUN(YP,Y,X,icase,IERR)  A user-supplied subroutine which evaluates
!   	the right hand side of the equation into vector YP.
!   	 The IERR is a parameter that may be set by FUN on a
!   	condition that requires the integration to terminate.
!
!
!
   REAL YW,YK1,YK2,YK3,YK4,YR1,YR2
   INTEGER  ND,icase
   DIMENSION YW(20),YK1(20),YK2(20),YK3(20),YK4(20),YR1(*),YR2(*)
        real yc(20),ycmag2,ydcmag2
!
!   		TYPE INITIAL VALUES.
!   WRITE(*,50)X0,(YR1(J),J=1,ND)
   DX2=DX/2.
!
!   		ADVANCE THE INTEGRATION.
   XHALF=X0+DX2
   XPLUS=X0+DX
   CALL RKFUN(YK1,YR1,X0,icase,IERR)
        ycmag2=0.
   do J=1,ND
           YK1(J)=YK1(J)*DX
           yc(j)=yk1(j)
           ycmag2=ycmag2+yc(j)**2
           YW(J)=YR1(J)+YK1(J)/2.
        enddo
   CALL RKFUN(YK2,YW,XHALF,icase,IERR)
   DO J=1,ND
           YK2(J)=YK2(J)*DX
           YW(J)=YR1(J)+YK2(J)/2.
        enddo
   CALL RKFUN(YK3,YW,XHALF,icase,IERR)
   DO J=1,ND
           YK3(J)=YK3(J)*DX
           YW(J)=YR1(J)+YK3(J)
        enddo
   CALL RKFUN(YK4,YW,XPLUS,icase,IERR)
        ydcmag2=0.
   DO J=1,ND
           YK4(J)=YK4(J)*DX
           ydcmag2=ydcmag2+(yk4(j)-yc(j))**2
           YR2(J)=YR1(J)+(YK1(J)+2.*YK2(J)+2.*YK3(J)+YK4(J))/6.
        enddo
        cond=ydcmag2/ycmag2
!
!   WRITE(*,50)X0,(Y(J),J=1,ND)
!50   FORMAT(8G10.3)
!
   RETURN
   END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Given a 2D function h(i,j) on a uniform grid, find its x and y gradients
! interpolated to position (x,y). Leading dimension of h is Lx+3>=nx+3.
! Gradients are piecewise linear. which causes a 2nd order magnitude
! reduction, but is fully adequate when h is a finite difference solution.
! Boundaries are specified by guard cells that go beyond allowed position.
! Positions are normalized such that x(0)=0, x(nx)=1. etc. So guard cells
! are x(-1)=-dx, x(nx+1)=1+dx, and allowed range is 0<=x<=1.
! Index will over-run if x is .ge. dx/2 past the allowed range.
! To use for a range xp0<xp<xpnx, pass x=(xp-xp0)/(xpnx-xp0) 
! and use returned gradient as dh/dxp=gx/(xpnx-xp0). (Et sim, y)

      subroutine grad2dint(h,Lx,nx,ny,x,y,gx,gy)
      integer Lx,nx,ny
      real h(-1:Lx+1,-1:ny+1),x,y,gx,gy

      im=int(x*nx)   ! Grid before
      ip=im+1        ! Grid after
      ig=nint(x*nx)  ! Nearest grid
      fx=x*nx-im     ! Transverse fraction
      fxg=x*nx-ig    ! Gradient fraction
!      write(*,*)ip,im,ig,x,nx

      jm=int(y*ny)
      jp=jm+1
      jg=nint(y*ny)
      fy=y*ny-jm
      fyg=y*ny-jg
!      write(*,*)jp,jm,jg,x,nx

      gxjp=(h(ig+1,jp)-h(ig,jp))*(0.5+fxg)                                  &
     &     +(h(ig,jp)-h(ig-1,jp))*(0.5-fxg)
      gxjm=(h(ig+1,jm)-h(ig,jm))*(0.5+fxg)                                  &
     &     +(h(ig,jm)-h(ig-1,jm))*(0.5-fxg)
      gx=(gxjm*(1.-fy)+gxjp*fy)*nx

      gyip=(h(ip,jg+1)-h(ip,jg))*(0.5+fyg)                                  &
     &     +(h(ip,jg)-h(ip,jg-1))*(0.5-fyg)
      gyim=(h(im,jg+1)-h(im,jg))*(0.5+fyg)                                  &
     &     +(h(im,jg)-h(im,jg-1))*(0.5-fyg)
      gy=(gyim*(1.-fx)+gyip*fx)*ny
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate function h(i,j) on uniform grid bilinearly. 
! See grad2dint above for domain scaling.
      real function bilin2d(h,Lx,nx,ny,x,y)
      integer Lx,nx,ny
      real h(-1:Lx+1,-1:ny+1),x,y

      im=int(x*nx)   ! Grid before
      ip=im+1        ! Grid after
      fx=x*nx-im     ! Transverse fraction
      jm=int(y*ny)
      jp=jm+1
      fy=y*ny-jm
      bilin2d=(h(im,jm)*(1-fx)+h(ip,jm)*fx)*(1-fy)                          &
     &     + (h(im,jp)*(1-fx)+h(ip,jp)*fx)*fy
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*********************************************************************


!***********************************************************************
subroutine orbitp
  use orbcartpoin
  character*120 string
  integer, dimension(2) :: imax
  real E(3)

  wpt=0.
  phip1=phiofrz(y0(1),0.)
  psi=phip1
  call getfield(y0,t,irktype,E) ! Initializes more conveniently.
! We enter this point with w0,Bsqpsi,r,th,z set, vth=0 by default
  if(iwritetype.eq.1)then
     Bsqpsi=Omega/sqrt(psi)
     B=Omega
!     write(*,*)'Bsqpsi',Bsqpsi,Omega,psi
  else
     B=Bsqpsi*sqrt(psi)
  endif
!  dt=0.047/B
  dt=dtxb/B
  if(lprint)write(*,'(a,f4.1,a,f5.2,a,f6.3,a,f5.2,a,f6.3,$)') &
       'w0=',w0,' Bsqpsi=',Bsqpsi,&
       ' psi=',psi,' r0=',y0(1),' B=',B !,' rho0=', y0(4)/B

  if(idoplots.ne.0)then
     if(lp)then
        call pfset(3)    ! Display and output plots
     else
        call pfset(-3)   ! Run continuously outputing plots
     endif
     call pltinit(-.99999,1.,-1.,0.)
     call charsize(.018,.018)
     call axis
     call axis2
     call axlabels('!Ac!@/!Ap!@','W!d!A|!@!d(r!dc0!d)/!Ay!@(r!dc0!d)')
  endif


  psic=psi
  do k=1,nwp
     t=0.
     if(nwp.le.2)then
        wp0=-(2*k-1)*wpf*psic
     else
!        wp0=-k*psi/(nwp+1.)
        wp0=-psic*(k/(nwp+1.))**aleph
     endif
! Initial position (zero velocity)
     y1=y0
! Initial r/theta velocities
     vmod0=sqrt(2.*(w0-wp0))
     y1(4)=cos(vangle0*3.1415926/180.)*vmod0
     y1(5)=sin(vangle0*3.1415926/180.)*vmod0
     pt0=y1(1)*y1(5)-B*y1(1)**2/2.   ! r*vt-B*r^2/2. = p_theta
! Approximate gyrocenter is one gyro-radius away from particle (y1)
! in frame rotating with ExB velocity v_tExB=-Er/Bz
     xg=y1(1)*cos(y1(2))-(y1(5)+E(1)/B)/B   ! r*cos(t) -(vt-vExB)/B
     yg=y1(1)*sin(y1(2))+y1(4)/B            ! r*sin(t) + vr/B
     rc0=sqrt(xg**2+yg**2)                            ! gyro-center radius
     rg0=sqrt(((y1(5)+E(1)/B)/B)**2+(y1(4)/B)**2)     ! gyro-radius
     rpmax=rc0+rg0
     rpmin=abs(rc0-rg0)
     rpave=(rpmax+rpmin)/2.
     phip1=phiofrz(rc0,y1(3))               ! at initial gc case.
!     phip1=phiofrz(rpave,y1(3))               ! at initial rpave case.

     if(k.eq.1)then
        psic=phiofrz(rc0,0.)
        if(lprint)write(*,'(a,f5.2,a,f5.3)')' rc0=',rc0,' psic=',psic
! Thermal rho=1/B.
        write(string,'(a,f5.3,a,f4.2,a,f5.3,a,f3.1,a,f4.2,a,f4.2)')&
             'B=',Omega,' W=',w0,' !Ay!@=',psi,' r!d0!d=',y0(1), &
             ' r!dc0!d=',rc0,' r!dg!d=',rg0
        
        if(idoplots.ne.0)call boxtitle(string(1:lentrim(string)))
        call winset(.true.)
     endif
     if(wp0+phip1.lt.0)then
        write(*,*)'NOT Starting with wp0',wp0,' below -phip1',-phip1
        write(*,*)'at rc0=',rc0,' z0=',y0(3)
        goto 12
     endif
     y1(6)=sqrt(2.*(wp0+phip1))

     if(lprint)write(*,'(i3,a,f8.4,a,f8.4,$)')k,' Wp0=',wp0,' vz=',y1(6)
     j=0
     wpp=wp0
     do i=1,nstep
        if(i.le.nplot.and.k.le.2)then
           r(i,k) =sqrt(y1(1)**2+y1(2)**2)
           th(i,k)=atan2(y1(2),y1(1))
           pt(i,k)=r(i,k)*(-y1(4)*sin(th(i,k))+y1(5)*cos(th(i,k))) &
                -B*r(i,k)**2/2.
           z(i,k) =y1(3)
           xplot(i,k)=r(i,k)*cos(th(i,k))
           yplot(i,k)=r(i,k)*sin(th(i,k))
           vr(i,k)=y1(4)
           vt(i,k)=y1(5)
           vz(i,k)=y1(6)
           phip(i,k)=phiofrz(r(i,k),y1(3))
           phir0(i,k)=phiofrz(rpave,y1(3))
!           Ep(i,k)=Eropsi*phip(i,k) Not correct when read in.
           w(i,k)=(vr(i,k)**2+vt(i,k)**2+vz(i,k)**2)/2.-phip(i,k)
           tv(i,k)=t
           imax(k)=i
        endif
   ! y1 etc are just (x,y,z,vx,vy,vz) and never change meaning.
        call RKADVC(dt,t,nvec,y1,irktype,y2,cond,IERR)
        phipi=phip1
        wpi=y2(6)**2/2.-phip1 ! Parallel energy at r0
        if(wpi*wpp.lt.0.)then
           if(lprint.and..false.)write(*,'(a,i6,a,f8.4,a,i4,a)') &
           ' Step=',i,' Wp sign flip' ,wpi,' after',j-1,' bounces'
           wpt=-wp0/psi
           wpp=wpi
!           exit
        endif
        if((.not.abs(y2(3)).lt.20.)  &
!             .or.(.not.abs(y2(1)).lt.2.*y0(1)) &  ! No good for low r.
             .or..not.wpi.lt.10.)then
           if(lprint)write(*,'(a,i6,a,f6.2,a,f8.4,a,i4,a)')&
           ' Step=',i,' z=',y2(3), &
           ' Wp=',wpi,' after',j-1,' bounces'
           wpt=-wp0/psi
           exit
        endif
        if(y2(3)*y1(3).le.0)then ! We crossed the z=0 center.
           if(j.eq.nbouncemax)exit
           j=j+1
           phip1=phiofrz(rc0,y2(3)) ! Evaluate new phi at initial gyrocenter
!           phip1=phiofrz(rpave,y2(3)) ! Evaluate new phi at rpave
           f1=abs(y2(3))/(abs(y1(3))+abs(y2(3)))
           f2=abs(y1(3))/(abs(y1(3))+abs(y2(3)))
           tc(j)=t

           wp(j)= f1*(y1(6)**2/2.-phipi) &  ! old
                 +f2*(y2(6)**2/2.-phip1)    ! new
           xi(j)=atan2(y1(5),y1(4)) ! This is absolute angle atan(vy/vx)
           xi(j)=xi(j)-atan2(y1(2),y1(1)) ! Subtract radial angle.
           if(wp(j).lt.-phip1.or.wp(j).gt.0..and.j.le.30)  &
            write(*,'(i5,a,f6.3,a,f8.4)')j,' xi=',xi(j),' wp/phip1=',wp(j)/phip1
        endif
        t=t+dt
        y1=y2
     enddo
     if((i.eq.nstep+1.or.j.eq.nbouncemax).and.lprint)then
        write(*,'(a,i4,$)')' j=',j
        write(*,73)' tb/2=',t/(j-1.), &
             ' tcyc=',2.*3.1415926/B,' ratio=',t/(j-1.)*B/2/3.1415926 
!             ,' rc0=',rc0
     endif
     if(idoplots.ne.0)call color(mod(k-1,14)+1)
     if(idoplots.ne.0)call polymark(xi/3.1415926,wp/phip1,j-2,3)
12   continue
  enddo
  if(idoplots.ne.0)call pltend
! This was the write used for the JGR plot
!  write(*,'(3f7.4,a)')Bsqpsi,wpt,Eropsi*sqrt(2.*(w0/psi-wpt)) &
!       ,'   b,wpt,Ervopsi3'
  if(iwritetype.eq.0)then
     write(*,'(4f10.4,a)')Bsqpsi,psi,wpt,y1(1),'  Bsqpsi,psi,wpt,rfinal'
  elseif(iwritetype.eq.1)then
     write(*,'(5f10.4,a)')Omega,psi,wpt,y1(1),Bsqpsi,'  Omega,psi,wpt,rfinal,Bsqp'
  endif
  do n=2,10,2
     wpr=((2.*B**2/n**2)**(-0.25)+1-2.**.25)**(-4)
  enddo
  
  if(nwp.le.2.and.idoplots.ge.2)call doplots2(imax,wp0,w0,pt0) 
  if(nwp.le.2.and.idoplots.ge.1)call doplots(imax,wp0,w0,pt0)  

!  583 format(i5,a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3)
  73 format(a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3)
end subroutine orbitp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine profilemin
  ! A modification of orbitmin so that the psi range is a range
  ! of radial positions. And we still find threshold Omega for confinement. 
  use orbcartpoin
  integer, parameter :: nbscan=10,nbiter=7
  real, dimension(npmmax) :: xioL,opp,ppp
  character*40 thestr,thewpf
  character*2 cel

  cel='L '
  OmegaLmax=Omega
  phip1=phiofrz(r0,0.)
  nwp=1
  iwritetype=2         ! No writing from orbitp
  idoplots=0

  call pfset(3)
  call pltinit(0.,yrfmax,0.,OmegaLmax)
  call charsize(.018,.018)
  call axis
  call axis2
  call axlabels('r',' ')
  call fwrite(wpm,iwidth,2,thewpf)
  call boxtitle('W!d!A|!@!d/!Ay!@='//thewpf(1:lentrim(thewpf)))
  do k=1,nerp
     wpf=1-(1-wpm)*k/nerp
     call fwrite(wpf,iwidth,2,thewpf)
     do i=1,npm   ! Here we have a range of r-positions.
        ri=(i-0.5)*yrfmax/npm
        psi=phiofrz(ri,0.)
        y0(1)=ri
        bmax=max(OmegaLmax,sqrt(psi))
!     write(*,*)ri,psi,wpf
2       do j=1,nbscan     ! Scan to find first confined orbit
           Omega=bmax*j/nbscan
           Bsqpsi=Omega/sqrt(psi)
           call orbitp
           if(wpt.eq.0)goto 1  ! A confined orbit.
        enddo
        bmax=bmax*1.5 !   If we did not find confined orbit, increment bmax
        goto 2
1       continue
        b2=bmax*(j)/nbscan
        b1=bmax*(j-1)/nbscan
        do j=1,nbiter          ! Iterate to refine the threshold omega.
           bb=(b2+b1)/2.
           Omega=bb
           Bsqpsi=Omega/sqrt(psi)
           call orbitp
           if(wpt.eq.0)then
              b2=bb
           else
              b1=bb
           endif
        enddo
        if(i.eq.1)rcr0=rc0
        write(*,'(i3,a,f8.5,a,f8.4,a,f8.4)')i,' psi=',psi,' Omega=',bb
        riplot(i,k)=ri
        psiplot(i,k)=psi
        omegaplot(i,k,1)=bb
     enddo
     call winset(.true.)
     call labeline(riplot,omegaplot(:,k,1),npm,thewpf,iwidth)
     call legendline(.7,.8,0,' !AW!@!dt!d')
     call dashset(4)
     call polyline([0.,rcr0],[1.,1.]*omegaplot(1,k,1),2)
     call drcstr(' r!dg0!d')
     call dashset(0)
  enddo
  call dashset(1)
  call color(1)
  call polyline(riplot,psiplot,npm)
  call legendline(.01,.03,0,' !Ay!@')
  call dashset(0)
  call pltend
end subroutine profilemin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program mainpoincare
  use orbcartpoin
  character*30 :: string
  real :: wn0=-10. ! Unset normalized w.
  do i=1,iargc()   ! Parse cmdline arguments.
     call getarg(i,string)
     if(string(1:1).ne.'-')read(string,*)filename
     if(string(1:2).eq.'-b')read(string(3:),*)Bsqpsi
     if(string(1:2).eq.'-O')then
        read(string(3:),*)Omega
        Bsqpsi=Omega/sqrt(psi)  ! Maybe not necessary now.
        iwritetype=1   ! Determines Omega setting has preference.
     endif
     if(string(1:2).eq.'-p')then
        if(Omega.ne.0.)then
           write(*,*)'DANGER: setting psi after Omega. Do NOT.',Omega
           stop
        endif
        read(string(3:),*)psi
     endif
     if(string(1:2).eq.'-q')then
        read(string(3:),*)psi
        psi=sqrt(psi)
     endif
     if(string(1:3).eq.'-nw')read(string(4:),*)nwp
     if(string(1:3).eq.'-ns')read(string(4:),*)nstep
     if(string(1:3).eq.'-mb')read(string(4:),*,end=103,err=103)wpm,psim,npm
     if(string(1:3).eq.'-mE')read(string(4:),*,end=103,err=103)nerp
103  continue
     if(string(1:2).eq.'-E')read(string(3:),*)Eropsi
     if(string(1:2).eq.'-r')read(string(3:),*)y0(1)
     if(string(1:2).eq.'-f')read(string(3:),*)wpf
     if(string(1:2).eq.'-W')read(string(3:),*)W0
     if(string(1:2).eq.'-w')read(string(3:),*)wn0
     if(string(1:2).eq.'-a')read(string(3:),*)aleph
     if(string(1:2).eq.'-v')read(string(3:),*)vangle0
     if(string(1:2).eq.'-c')lp=.not.lp
     if(string(1:2).eq.'-L')lo=.not.lo
     if(string(1:2).eq.'-i')lpi=.not.lpi
     if(string(1:2).eq.'-d')lprint=.not.lprint
     if(string(1:2).eq.'-h')goto 201
     if(string(1:2).eq.'-?')goto 201
  enddo
  if(nerp.gt.nerpmax)stop 'nerp specified too big'
  if(wn0.gt.-10.)w0=psi*wn0
  if(Omega.eq.0)Omega=Bsqpsi*sqrt(psi)

  if(wpm.eq.0)then
     call orbitp
  else            ! Find minimum confined
     call profilemin
  endif
  call exit
  
201 continue
  write(*,*)'Usage orbcartpoin [flags] [filename]'
  write(*,*)'-b... B/sqpsi, -p... psi, -q... sqrt(psi), -E... Er/psi, -r... r0'
  write(*,*)'-W... W0 (total energy), -w ... w0 (normalized alternate)'
  write(*,*)'-O... Omega, resets actual B, B/sqpsi value; use after all -p'
  write(*,*)'-nw... N-energies [if 1 plot orbits], -f... the single energy'
  write(*,*)'-a... set aleph power to concentrate wp values near zero. '
  write(*,*)'-v... set vangle0 degrees of initial velocity to direction r/y.'
  write(*,*)'-ns... N-steps, -c run continuously, -d do not write, -h help'

  write(*,*)'-mb[<wpm>,<psim>,<npm>] Find the minimum B that confines the orbit'
  write(*,*)'      of fractional depth wpm, up to psimax [0.5], in npm steps [100]'
  write(*,*)'-mE   set number of different cases (Eropsi,wpf) for -mb call.'
  write(*,*)'-L    toggle plotting disruption/resonance scaling'
  write(*,*)'-i    toggle plotting input phi contours'
end program mainpoincare

