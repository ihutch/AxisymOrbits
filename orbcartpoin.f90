  ! Integrate orbits in a specified potential to generate a poicare plot
  ! Version that calculates the orbit in cartesian coordinates.
  ! And uses a file read in to specify the potential and E-field.
module orbcartpoin
  real :: Bsqpsi=.85, Eropsi=.01, psi=1., Omega=0.
  real :: r0,rc0,rg0,rpmin,rpmax,rpave,wpf=0.5,B,w0=1.,aleph=1.,wpt
  real :: wpm=0.,psim=0.5   ! Minimum confined energy, for psimax
  real :: dtxb=0.04363     ! dt x Omega=2pi/(steps per cyclotron period=144)
  real :: yrfmax,zrfmax,psiave,tk
  integer, parameter :: npmmax=100,nerpmax=20
  real, dimension(npmmax,nerpmax) :: psiplot,riplot,wpplot,wmeanplot,psimeanplot
  real, dimension(npmmax,nerpmax,2) :: omegaplot
  integer :: iotype=0, nerp=1,npm=11  ! type, n Eropsi's, psi divisions
  logical :: lp=.true.,lprint=.true.,lo=.true.,lpi=.false.
  logical :: lgread=.false.,lhread=.false.,lpg=.true.
  integer, parameter :: nbouncemax=200  ! Max bounces to store.
  integer, parameter :: nplot=25000     ! Max steps to store.
  integer, parameter :: nvec=6
  integer :: iwritetype=0,ibtype=0,idoplots=2,irktype=0,idebug=0
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
      real :: zrf,yrf,zs,rh,zsa,gy,gz
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
      if(ys.lt.0)write(*,*)'phiofrz rh negative',rh
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
      call minmax(r(:,np),imax,rmin1,rmax)
      rmin=rmin1-.3
      rmax=rmax+.3
      call minmax(vz(:,np)**2/2.-phip(:,np),imax(np),wmin,wmax)

      nbb=nint(3.1415926/dtxb)
      rmean=0
      wmean=0
      ncyco=2*3
      do i=1,ncyco*nbb
         rmean=rmean+r(i,1)
         wmean=wmean+vz(i,1)**2/2.-phip(i,1)
      enddo
      rmean=rmean/(ncyco*nbb)
      wmean=wmean/(ncyco*nbb)
      write(*,*)'rc0,rmean,Wpmean',rc0,rmean,wmean

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
      call color(2)
      call polyline(xptemp(nbb),yptemp(nbb),imax-2*nbb)
      call winset(.false.)
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
         psicplot(i)=-phiofrz(rmin1,zlplot(i))    ! At lowest r.
!         psilplot(i)=-phiofrz(r0,zlplot(i))      ! At initial r.
         psilplot(i)=-phiofrz(rc0,zlplot(i))      ! At initial gyrocenter.
      enddo
      call pltinit(zmin,zmax,psicplot(nzplot/2), &
           max(wmax,-psicplot(nzplot/2)*.2))
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
! Set B if necessary from phip1, and write out parameters.
  call getfield(y0,t,irktype,E) ! Field at initial point.
  phip1=phiofrz(y0(1),0.)
  if(ibtype.ne.1)then
     B=Bsqpsi*sqrt(phip1)   ! Bsqpsi at point determines B.
  else
     Bsqpsi=B/sqrt(phip1)   ! B at point determines Bsqpsi.
  endif
  dt=dtxb/B
  if(idoplots.ne.0)then  ! Set up Poincare plot
     if(lp)then
        call pfset(3)    ! Display and output plots
     else
        call pfset(-3)   ! Run continuously outputing plots
     endif
     call pltinit(-.99999,1.,-1.,max(0.,-2.*wpf))
     call charsize(.018,.018)
     call axis
     call axis2
     call axlabels('!Ac!@/!Ap!@','W!d!A|!@!d(r!dc0!d)/!Ay!@(r!dc0!d)')
  endif


  psic=phip1 ! Will be updated in a moment.
  do k=1,nwp
     t=0.
! Set parallel energy
     if(nwp.le.2)then
        wp0=-(2*k-1)*wpf*psic
     else
        wp0=-psic*(k/(nwp+1.))**aleph
     endif
! Initial position (zero velocity)
     y1=y0
! Initial r/theta velocities
     vmod0=sqrt(2.*(w0-wp0))            ! Should add potential? No.
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
        psic=phiofrz(rc0,0.)    ! Determines Wp range.
        if(lprint)write(*,'(a,f4.1,a,f5.2,a,f5.3,a,f6.3,a,f6.3,a,f5.2,a,f5.3)')&
             'w0=',w0,' Bsqpsi=',Bsqpsi,' r0=',y0(1),&
             ' psic=',psic,' B=',B,' rc0=',rc0
        write(string,'(a,f5.3,a,f4.2,a,f5.3,a,f3.1,a,f4.2,a,f4.2)')&
             'B=',Omega,' W=',w0,' !Ay!@=',psic,' r!d0!d=',y0(1), &
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
           w(i,k)=(vr(i,k)**2+vt(i,k)**2+vz(i,k)**2)/2.-phip(i,k)
           tv(i,k)=t
           imax(k)=i
        endif
   ! y1,2 are just (x,y,z,vx,vy,vz) and never change meaning.
        call RKADVC(dt,t,nvec,y1,irktype,y2,cond,IERR)
        wpi=y2(6)**2/2.-phiofrz(rc0,y1(3)) ! Parallel energy at rc0
        if((.not.abs(y2(3)).lt.20.)  & ! Untrapped.
             .or..not.wpi.lt.10.)then
           if(lprint)write(*,'(a,i6,a,f6.2,a,f8.4,a,i4,a)')&
           ' Step=',i,' z=',y2(3), &
           ' Wp=',wpi,' after',j-1,' bounces'
           wpt=-wp0/psic
           exit
        endif
        if(y2(3)*y1(3).le.0)then ! Crossed the z=0 center. Document.
           if(j.eq.nbouncemax)exit
           j=j+1
           phipi=phip1
           phip1=phiofrz(rc0,y2(3)) ! Evaluate new phi at initial gyrocenter
!           phip1=phiofrz(rpave,y2(3)) ! Evaluate new phi at rpave
           f1=abs(y2(3))/(abs(y1(3))+abs(y2(3)))
           f2=abs(y1(3))/(abs(y1(3))+abs(y2(3)))
           tc(j)=t

           wp(j)= f1*(y1(6)**2/2.-phipi) &  ! old
                 +f2*(y2(6)**2/2.-phip1)    ! new
           xi(j)=atan2(y1(5),y1(4)) ! This is absolute angle atan(vy/vx)
           xi(j)=xi(j)-atan2(y1(2),y1(1)) ! Subtract radial angle.
!           if(wp(j).lt.-phip1.or.wp(j).gt.0..and.j.le.30)  &
!           write(*,'(i5,a,f6.3,a,f8.4)')j,' xi=',xi(j),' wp/phip1=',wp(j)/phip1
        endif
        t=t+dt
        y1=y2
     enddo
     if((i.eq.nstep+1.or.j.eq.nbouncemax).and.lprint)then ! Final print.
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
  if(iwritetype.eq.0)then
     write(*,'(4f10.4,a)')Bsqpsi,psic,wpt,y1(1),'  Bsqpsi,psi,wpt,rfinal'
  elseif(iwritetype.eq.1)then
     write(*,'(5f10.4,a)')Omega,psic,wpt,y1(1),Bsqpsi,'  Omega,psi,wpt,rfinal,Bsqp'
  endif
  
  if(nwp.le.2.and.idoplots.ge.2)call doplots2(imax,wp0,w0,pt0) 
  if(nwp.le.2.and.idoplots.ge.1)call doplots(imax,wp0,w0,pt0)  

  73 format(a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3)
end subroutine orbitp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine orbitc
! Specifies gyrocenter rather than initial position, and takes the
! vangle and W0 to be in the drift-rotating frame using E(rc)/B
! In this routine, psic is the potential at rc0,0; psiave is the 
! potential second-order gyro averaged (cartesian) about rc0,0;
! and phiave is the gyroaveraged potential about rc0,z (close to z=0)
  use orbcartpoin
  character*120 string
  integer, dimension(2) :: imax
  real E(3)
  psic=0.
  wpt=0.
! Set E and B if necessary from 
  call getfield(y0,t,irktype,E) ! Field at initial point.
  psic=phiofrz(y0(1),0.)  ! y0(1) is gyrocenter, can't average yet.
  psiave=psic
  do i=1,2
     if(ibtype.ne.1)then
        B=Bsqpsi*sqrt(psiave)   ! Bsqpsi at point determines B.
     else
        Bsqpsi=B/sqrt(psiave)   ! B at point determines Bsqpsi.
     endif
! Orbit advance uses only B, not Bsqpsi.
! Semi iterative 2nd order estimate of gyroaveraged phi, at rc0
     rg0=sqrt(2.*W0)/B ! Not including wp0 or psiave.
     psiave=(2*phiofrz(y0(1),0.) &
          +phiofrz(y0(1)+rg0,0.)+phiofrz(abs(y0(1)-rg0),0.))/4.
!     write(*,*)'psic,psiave',psic,psiave,B,Bsqpsi
  enddo
  phiave=psiave ! Barely necessary.
  
  dt=dtxb/B
  if(idoplots.ne.0)then  ! Set up Poincare plot
     if(lp)then
        call pfset(3)    ! Display and output plots
     else
        call pfset(-3)   ! Run continuously outputing plots
     endif
     call pltinit(-.99999,1.,-1.,max(0.,-2.*wpf))
     call charsize(.018,.018)
     call axis
     call axis2
     call axlabels('!Ac!@/!Ap!@','W!d!A|!@!d(r!dc0!d)/!Ay!@(r!dc0!d)')
  endif

  do k=1,nwp
     tk=0.
! Set parallel energy
     if(nwp.le.2)then
        wp0=-(2*k-1)*wpf*psiave
     else
        wp0=-psic*(k/(nwp+1.))**aleph
     endif
! Initial position
     y1=y0   
! Initial x/y velocities in the drift frame.
!     vmod0=sqrt(2.*(psic+w0-wp0))  ! Should add potential? Not if wp0 includes.
!     vmod0=sqrt(2.*(w0-wp0))  ! Should add potential? No. wp0 includes.
! The above were my first tries and neither are consistent.

! Here wp0 is parallel kinetic energy minus psiave.  The perpendicular
! kinetic energy in the drift frame must be set to the total energy in the
! fixed frame W0 minus the parallel kinetic energy, plus the actual
! potential at point, corrected for the drift energy.  But if the drift 
! energy is not large, perhaps it can be ignored. If so,
! vmod0^2/2=W0-(wp0+psiave)+psipoint. However point is not known
! because the position depends on rg0 which depends on vmod. So iterate
     psip=psic
     do i=1,3
        vmod0=sqrt(2.*(W0+psip-(wp0+psiave)))
        y1(4)=cos(vangle0*3.1415926/180.)*vmod0
        y1(5)=sin(vangle0*3.1415926/180.)*vmod0
        rg0=vmod0/B                       ! gyro-radius
! If gyrocenter is at x=y1(1), y=0, then
        xg=y1(1)
        yg=0.
        rc0=xg                            ! gyro-center radius
! and the particle is actually at
        y1(1)=xg+y1(5)/B
        y1(2)=yg-y1(4)/B
        psip=phiofrz(sqrt(y1(1)**2+y1(2)**2),0.)
!        if(lprint)write(*,*)rg0,xg,psip,(y1(4)**2+y1(5)**2)/2.
     enddo
     rpmax=rc0+rg0
     rpmin=abs(rc0-rg0)
     rpave=(rpmax+rpmin)/2.                         ! =rc0 here.
! Transform initial velocities to rest frame by adding drift (at gc)
     y1(4)=y1(4)+E(2)/B  ! Should not be needed because E(2)=0
     y1(5)=y1(5)-E(1)/B
! Now in rest-frame
    ! r*vt-B*r^2/2. = p_theta
     pt0=y1(1)*y1(5)-y1(2)*y1(4)-B*(y1(1)**2+y1(2)**2)/2.
     psic=(2.*phiofrz(rc0,y1(3))+phiofrz(rpmax,y1(3))+phiofrz(rpmin,y1(3)))/4.
     
     if(k.eq.1)then
        if(lprint)write(*,'(a,f5.2,a,f5.2,a,f6.3,a,f6.3,a,f6.3,a,f5.2,a,f5.2)')&
             'w0=',w0,' Bsqpsi=',Bsqpsi,' x0=',y1(1),&
             ' psic=',psic,' B=',B,' rc0=',rc0,' rg0=',rg0
        write(string,'(a,f5.3,a,f4.2,a,f5.3,a,f4.1,a,f4.2,a,f4.2)')&
             'B=',Omega,' W=',w0,' !Ay!@=',psic, & !' x!d0!d=',y1(1), &
             ' r!dc0!d=',rc0,' r!dg!d=',rg0
        
        if(idoplots.ne.0)call boxtitle(string(1:lentrim(string)))
        call winset(.true.)
     endif
     if(wp0+psiave.lt.0)then
        write(*,*)'NOT Starting with wp0',wp0,' below -psiave',-psiave
        write(*,*)'at rc0=',rc0,' z0=',y0(3),' wpf=',wpf
        goto 12
     endif
     y1(6)=sqrt(2.*(wp0+psiave))
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
           phir0(i,k)=phiofrz(rpave,y1(3)) ! Ought to be gyroaverage.
           w(i,k)=(vr(i,k)**2+vt(i,k)**2+vz(i,k)**2)/2.-phip(i,k)
           tv(i,k)=tk
           imax(k)=i
        endif
   ! y1,2 are just (x,y,z,vx,vy,vz) and never change meaning.
        call RKADVC(dt,t,nvec,y1,irktype,y2,cond,IERR)
        wpi=y2(6)**2/2.-phiofrz(rc0,y1(3)) ! Parallel energy at rc0
        if((.not.abs(y2(3)).lt.20.)  & ! Untrapped.
             .or..not.wpi.lt.10.)then
           if(lprint)write(*,'(a,i6,a,f6.2,a,f8.4,a,i4,a)')&
           ' Step=',i,' z=',y2(3), &
           ' Wp=',wpi,' after',j-1,' bounces'
           wpt=-wp0/psic
           exit
        endif
        if(y2(3)*y1(3).le.0)then ! Crossed the z=0 center. Document.
           if(j.eq.nbouncemax)exit
           j=j+1
           phiavei=phiave
           phiave=(2.*phiofrz(rc0,y2(3))+phiofrz(rpmax,y2(3)) &
                +phiofrz(rpmin,y2(3)))/4.  ! Average phi.
           f1=abs(y2(3))/(abs(y1(3))+abs(y2(3)))
           f2=abs(y1(3))/(abs(y1(3))+abs(y2(3)))
           tc(j)=t

           wp(j)= f1*(y1(6)**2/2.-phiavei) &  ! old
                 +f2*(y2(6)**2/2.-phiave)    ! new
           xi(j)=atan2(y1(5),y1(4)) ! This is absolute angle atan(vy/vx)
           xi(j)=xi(j)-atan2(y1(2),y1(1)) ! Subtract radial angle.
        endif
        tk=tk+dt
        y1=y2
     enddo
     if((i.eq.nstep+1.or.j.eq.nbouncemax).and.lprint)then ! Final print.
        write(*,'(a,i4,$)')' j=',j
        write(*,73)' tb/2=',t/(j-1.), &
             ' tcyc=',2.*3.1415926/B,' ratio=',t/(j-1.)*B/2/3.1415926 
!             ,' rc0=',rc0
     endif
     if(idoplots.ne.0)call color(mod(k-1,14)+1)
     if(idoplots.ne.0)call polymark(xi/3.1415926,wp/psiave,j-2,3)
12   continue
  enddo
  if(idoplots.ne.0)call pltend
  if(iwritetype.eq.0)then
     write(*,'(4f10.4,a)')Bsqpsi,psic,wpt,y1(1),'  Bsqpsi,psi,wpt,rfinal'
  elseif(iwritetype.eq.1)then
     write(*,'(5f10.4,a)')Omega,psic,wpt,y1(1),Bsqpsi,'  Omega,psi,wpt,rfinal,Bsqp'
  endif
  
  if(nwp.le.2.and.idoplots.ge.2)call doplots2(imax,wp0,w0,pt0) 
  if(nwp.le.2.and.idoplots.ge.1)call doplots(imax,wp0,w0,pt0)  

  73 format(a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3)

end subroutine orbitc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine profilemin
  ! A modification of orbitmin so that the psi range is a range
  ! of radial positions. And we still find threshold Omega for confinement. 
  use orbcartpoin
  integer, parameter :: nbscan=10,nbiter=7
!  real, dimension(npmmax) :: xioL
  character*40 thewpf
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
!           write(*,*)'orbit call',j,Omega,Bsqpsi,lpg
           if(lpg)call orbitc
           if(.not.lpg)call orbitp
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
           if(lpg)call orbitc
           if(.not.lpg)call orbitp
           if(wpt.eq.0)then
              b2=bb
           else
              b1=bb
           endif
        enddo
        if(i.eq.1)rcr0=rc0
        write(*,'(i3,a,f8.4,a,f8.5,a,f8.4,a,f8.4)') &
             i,' rc0=',rc0,' psi=',psi,' Omega=',bb
        riplot(i,k)=ri
        psiplot(i,k)=psi
        omegaplot(i,k,1)=bb
     enddo
     call winset(.true.)
     call labeline(riplot,omegaplot(:,k,1),npm,thewpf,iwidth)
     call legendline(.7,.8,0,' !AW!@!dt!d')
  enddo
  call dashset(1)
  call color(1)
  call polyline(riplot,psiplot,npm)
  call legendline(.01,.03,0,' !Ay!@')
  call dashset(0)
  call pltend
end subroutine profilemin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wploss
  use orbcartpoin
! Search for the parallel energy loss limit profile over radius.
  integer :: nwiter=5,nOm
  real, dimension(npmmax) :: rpmplot,phipmplot
  real Omx
  character*20 string
  nwp=1                ! Single orbits
  iwritetype=2         ! No writing from orbitc
  idoplots=0           ! No plotting from orbitc
  phip1=phiofrz(r0,0.) ! Initialize potential and grid arrays.
  call charsize(0.018,0.018)
  if(wpm.le.1)then
     write(*,*)'wploss using wpm=',1.1
     wpm=1.1
  endif
  
  do i=1,npmmax        ! Radial potential profile
     rpmplot(i)=yrfmax*(i-1)/(npmmax-1)
     phipmplot(i)=phiofrz(rpmplot(i),0.)
  enddo
  Omx=Omega
  nOm=nerp
  do k=1,nOm                  ! For range of Omega
     Omega=omx*(k+1)/(nOm+1)  ! Bigger gap at bottom.
     B=Omega
     nradius=11
!     write(*,*)'nOm,k,Omega,Omx',nOm,k,Omega,Omx
     do i=1,nradius            ! Radial gyrocenter position
        y0(1)=yrfmax*(i-1)/(nradius-1)
        psi=phiofrz(y0(1),0.)
        if(ibtype.ne.1)then; write(*,*)'You need switch -O...'; stop;
        endif
        lprint=.true.
        do j=1,npm             ! wpsearch upward
           wpf=1.-wpm*j/npm
!           write(*,*)'npm,wpm,j,wpf=',npm,wpm,j,wpf
           if(j.gt.1)lprint=.false.
           call orbitc         ! Follow the wpf orbit.
           if(wpt.ne.0)goto 1  ! An unconfined orbit.
        enddo
        write(*,*)'Did not find unconfined orbit up to',-wpf
        exit
1       continue
        ilowun=0
        if(j.eq.1)then
           write(*,'(a,3i3,6f8.4)')'Lowest orbit is unconfined',&
                i,j,npm,wpf,wpm,psiave,r(1,1),vz(1,1)
           ilowun=1
        else
     ! Now j is unconfined, and j-1 is confined.
           w2=wpf
           w1=1-(1-wpm)*(j-1)/npm
           do j=1,nwiter       ! Iterative solution refinement.
              wpf=(w1+w2)/2.
              call orbitc
              if(wpt.eq.0)then; w1=wpf; else; w2=wpf
              endif
           enddo
           if(wpt.ne.0)then; wpf=w1; call orbitc ! Last always confined.
           endif
        endif
! Now we have just run the highest confined orbit or have found that the
! lowest is unconfined. tk is now the last time before the end of orbit.
        nbb=nint(3.1415926/dtxb)  ! Number of steps in half a gyroperiod
! We should take averages over some number of half-gyroperiods low enough
! not to go past the end of an unconfined lowest orbit. This limit is
!       ngby2*nbb*dt=tk i.e. ngby2=tk/(3.1415926/B)
        ngby2=int(tk/(3.1415926/B))
        ncyco=min(10,ngby2)
        navsteps=ncyco*nbb
        rmean=0
        psimean=0
        wmean=0
        do ic=1,navsteps
           rmean=rmean+r(ic,1)
           psimean=psimean+phiofrz(r(ic,1),0.)
           wmean=wmean+(vz(ic,1)**2/2.-phip(ic,1))
        enddo
        rmean=rmean/(ncyco*nbb)
        psimean=psimean/(ncyco*nbb)
        wmean=wmean/(ncyco*nbb)/psimean
        if(ilowun.eq.1)wmean=-wpf  ! Not a confined orbit; plot low.
        write(*,'(2i3,a,f6.3,a,f6.3,a,f6.3,a,f8.4,a,f8.4)') &
             k,i,' rc0=',rc0,' rg0=',rg0,' rmean=',rmean,&
             ' wmean=',wmean,' wpf=',wpf
        riplot(i,k)=y0(1)
        psiplot(i,k)=phiofrz(y0(1),0.)
        wpplot(i,k)=wpf
        wmeanplot(i,k)=wmean
        psimeanplot(i,k)=psimean
!        write(*,*)'psiave,psimean',psiave,psimean
     enddo
     if(k.eq.1)then
        wpmin=-0.6
        call pltinit(0,yrfmax,wpmin,.1)
        call axis
        call axlabels('r!dc!d','!p!o-!o!qW!d!A|!@!d/!p!o-!o!q!Ay!@')
        call winset(.true.)
     endif
     call color(k)
     call fwrite(B/sqrt(phiofrz(0.,0.)),iwidth,1,string)
     call dashset(k)
     call polyline(riplot(:,k),wmeanplot(:,k),i-1)
     call polyline([0.,rg0],(1-[1.,1.]*k/nOm*.3)*wpmin,2)
     call drcstr(' '//string(1:iwidth))
     call accisflush
  enddo   ! End of Omega iteration
  call color(15)
  call dashset(0)
  call winset(.false.)
!  call jdrwstr(wx2nx(0.05),wy2ny(0.06),'Labels: !AW!@/!Ay!@!u1/2!u',1.)
  call jdrwstr(wx2nx(rg0),wy2ny(0.63*wpmin),'!AW!@/!Ay!@!d0!d!u1/2!u',1.)
  call jdrwstr(wx2nx(0.),wy2ny(0.12+wpmin),'r!dgt!d:',-1.2)
  call axptset(1.,0)
  call scalewn(0.,yrfmax,0.,1.3*phiofrz(0.,0.),.false.,.false.)
  call ticrev
  call altyaxis(1.,1.)
  call axptset(0.,0.)
  call ticrev
  call polyline(rpmplot,phipmplot,npmmax)
  ix=int(0.6*npmmax)
  call jdrwstr(wx2nx(rpmplot(ix)),wy2ny(phipmplot(ix))+.02, &
       '!p!u-!u!q!Ay!@(r)!A_!@',1.)
!  call pltend
! psimeanplot additional?
!  call pltinit(0.,yrfmax,0.,2*phiofrz(0.,0.))
!  call axis; call axis2
!  call polyline(rpmplot,phipmplot,npmmax)
!  ix=int(0.8*npmmax)
!  call jdrwstr(wx2nx(rpmplot(ix)),wy2ny(phipmplot(ix))+.02,'!Ay!@(r)!A_!@',1.)
  do k=1,nOm
!     call color(k)
     call dashset(k)
     call polyline(riplot(:,k),psimeanplot(:,k),nradius)
  enddo
  
  call pltend
  
end subroutine wploss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        Bsqpsi=Omega/sqrt(psi)
        B=Omega
        iwritetype=1   ! Determines Omega setting has preference.
        ibtype=1
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
     if(string(1:3).eq.'-mb')then
        iotype=1
        read(string(4:),*,end=103,err=103)wpm,psim,npm
     endif
     if(string(1:3).eq.'-mw')then
        iotype=2
        read(string(4:),*,end=103,err=103)wpm,psim,npm
     endif
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
     if(string(1:2).eq.'-g')lpg=.not.lpg
     if(string(1:2).eq.'-d')lprint=.not.lprint
     if(string(1:2).eq.'-h')goto 201
     if(string(1:2).eq.'-?')goto 201
  enddo
  if(nerp.gt.nerpmax)stop 'nerp specified too big'
  if(wn0.gt.-10.)w0=psi*wn0
  if(Omega.eq.0)Omega=Bsqpsi*sqrt(psi)
  if(lp)then; call pfset(3); else; call pfset(-3);
  endif
  if(iotype.eq.0)then
     if(lpg)call orbitc
     if(.not.lpg)call orbitp
  elseif(iotype.eq.1)then            ! Find minimum confined
     call profilemin
  elseif(iotype.eq.2)then
     call wploss
  endif
  call exit
  
201 continue
  write(*,*)'       Usage orbcartpoin [-flag -flag ...] [filename]'
  write(*,6)'Flags, showing current value: ['
  write(*,7)'-b... B/sqpsi, -p... psi, -q... sqrt(psi), -E... Er/psi, -r... r0['
  write(*,'(5f13.3)')Bsqpsi,psi,sqrt(psi),Eropsi,y0(1)
  write(*,7)'-W... W0 (total energy), -w ... w0 (normalized alternate)[',W0
  write(*,7)'-O... Omega, resets actual B, B/sqpsi value; use after -p[',Omega
  write(*,8)'-nw... N-energies [if =1, plot orbits]                   [',nwp
  write(*,7)'-f...  the single energy when nw=1                       [',wpf
  write(*,7)'-a... set aleph power to concentrate wp values near zero.[',aleph
  write(*,7)'-v... set vangle0 degrees of velocity to direction r/y.  [',vangle0
  write(*,'(a,i8)')'-ns... N-steps [',nstep
  write(*,7)'-mb[<wpm>,<psim>,<npm>] Find the minimum B that confines the orbit'
  write(*,8)'   of fractional depth wpm, up to psimax, in npm steps  [',npm
  write(*,8)'-mE... set number of cases (Eropsi or B) for -mb/w call.[',nerp
  write(*,7)'-mw[<wpm>,<psim>,<npm>] Find max wpm that confines      [',wpm,psim
  write(*,9)'-L    toggle plotting disruption/resonance scaling      [',lo
  write(*,9)'-i    toggle plotting input phi contours                [',lpi
  write(*,9)'-g    toggle gyrocenter radius initialization.          [',lpg
  write(*,6)'filename non-switch names potential to read             [ ', &
       filename(1:lentrim(filename))
  write(*,7)'-c run continuously, -d do not write, -h help'
6 format(6a)  
7 format(a,2f7.3)
8 format(a,i4)
9 format(a,l4)
end program mainpoincare

