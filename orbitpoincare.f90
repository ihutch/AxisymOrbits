  ! Integrate orbits in a specified potential to generate a poicare plot
module orbitpoincare
  real :: Bsqpsi=.85, Eropsi=.01, psi=1., Omega=0.
  real :: r0,wpf=0.5,B,w0=1.,alpha=1.,wpt
  real :: wpm=0.,psim=0.5   ! Minimum confined energy, for psimax
  integer, parameter :: npmmax=100,nerpmax=20
  real, dimension(npmmax,nerpmax) :: psiplot
  real, dimension(npmmax,nerpmax,2) :: omegaplot
  integer :: nerp=1,npm=npmmax   ! n Eropsi's, psi divisions of psimax
  logical :: lp=.true.,lprint=.true.,lo=.true.
  integer, parameter :: nbouncemax=200  ! Max bounces to store.
  integer, parameter :: nplot=25000     ! Max steps to store.
  integer, parameter :: nvec=6
  integer :: iwritetype=0,idoplots=2
  integer :: nwp=39,nstep=1000000       ! Default maximum steps to take.
  real, dimension(nplot,2) :: r,th,z,vr,vt,vz,w,tv,pt,xplot,yplot,phip,Ep
  real, dimension(nbouncemax) :: wp,xi,tc
  real, dimension(nvec) :: y0,y1,y2
  real ::  vangle0=0.          ! Angle of initial v to r-direction degrees.
  data y0/10.,0.,0.,0.,0.,0./ !Initial orbit (position) defaults
contains
!***********************************************************************
   subroutine RKFUN(yp,y,t,icase,ierr)
     real yp(6),y(6),t
     integer icase,ierr
! Lorentz RK advance function. Evaluates derivative of y into yp.
! Cylindrical coordinates r theta z, vr, vt, vz
     real E(3)  ! Er,Et(=0),Ez
! Evaluate E.
     call getfield(y,t,icase,E)
! Enter derivatives.
      yp(1)=y(4)        ! vr
      yp(2)=y(5)/y(1)   ! vt/r = d/dt theta
      yp(3)=y(6)        ! vz
      yp(4)=y(5)**2/y(1)-E(1)-y(5)*B   ! vt^2/r -Er-vt*B
      yp(5)=-y(4)*y(5)/y(1)+y(4)*B    ! -vr*vt/r+vr*B 
      yp(6)=-E(3)
      ierr=0
    end subroutine RKFUN
!***********************************************************************
    subroutine getfield(y,t,icase,E)
      real y(6),t,E(3)
      integer icase
      psir=psi*max((1.+(r0-y(1))*Eropsi),0.) ! Linear 
!      psir=psi*exp((r0-y(1))*Eropsi)         ! Exponential psi(r) 
      sc=1./cosh(y(3)/4.)
      E(1)=0.
      if(psir.gt.0.)E(1)=Eropsi*psi*sc**4    ! Linear
!      E(1)=Eropsi*psir*sc**4                 !Exponential
      E(2)=0.
      E(3)=psir*sinh(y(3)/4.)*sc**5
    end subroutine getfield
!***********************************************************************
    real function psiofrz(rh,zh)
      cc=1./cosh(zh/4.)
      psir=psi*max((1.+(r0-rh)*Eropsi),0.)      ! Linear
!      psir=psi*exp((r0-rh)*Eropsi)              ! Exponential
      psiofrz=cc**4*psir
    end function psiofrz
!***********************************************************************
    subroutine doplots2(imax,wp0,w0,pt0)
      integer, dimension(2) :: imax
      integer, parameter :: nzplot=100
      real, dimension(nzplot) :: zplot,psiplot
      if(lp)then; call pfset(3); else; call pfset(-3);
      endif

      np=min(nwp,2)
      call multiframe(3,1,0)
      zmin=-10.
      zmax=10.
      call minmax(r(:,np),imax,rmin,rmax)
      rmin=rmin-.3
      rmax=rmax+.3
      call minmax(vz(:,np)**2/2.-phip(:,np),imax(np),wmin,wmax)

      call pltinit(zmin,zmax,rmin,rmax)
      call ticnumset(4)
      call charsize(.015,.015)
      call axis
      call jdrwstr(.08,wy2ny(r0+.5),'r',-1.)
      call winset(.true.)
      call color(4)
      call polymark(z,r,1,1)
      call polyline(z,r,imax)
      call color(15)
!      call pltend

!      call autoplot(z,vz**2/2.-phip,imax)
      call pltinit(zmin,zmax,wmin,0.02)
      call charsize(.015,.015)
      call axis
      call jdrwstr(.08,wy2ny(-0.6*wpf),'W!d!A|!@!d',-1.)
      call winset(.true.)
      call color(4); call polyline(z,vz**2/2.-phip,imax)
      call polymark(z,vz**2/2.-phip,1,1)
      call color(1);
!      write(*,*)'imax=',imax
      call polyline(z(:,2),vz(:,2)**2/2.-phip(:,2),imax)
      call polymark(z(:,2),vz(:,2)**2/2.-phip(:,2),1,1)
      call color(15)
      call dashset(4)
      call polyline([zmin,zmax],[0.,0.],2)
      call dashset(0)

      do i=1,nzplot
         zplot(i)=zmin+(zmax-zmin)*(i-1.)/(nzplot-1)
         psiplot(i)=-psi/cosh(zplot(i)/4.)**4
      enddo
      call pltinit(zmin,zmax,-psi-.1,0.1)
      call charsize(.015,.015)
      call axis
      call axlabels('z','')
      call axptset(1.,0.)
      call ticrev; call altyaxis(Eropsi,Eropsi); call ticrev
      call axptset(0.,0.)
      call jdrwstr(.08,wy2ny(-.2*psi),'W!d!A|!@!d',-1.)
      call jdrwstr(.12,wy2ny(-.2*psi),'-!Af!@',1.)
      call jdrwstr(.96,wy2ny(-.2*psi),'-E!dr!d',-1.)
      call winset(.true.)
      call polyline(zplot,psiplot,nzplot)
      call color(1);
      call polyline(z(:,2),vz(:,2)**2/2.-phip(:,2),imax)
!      call polymark(z(:,2),vz(:,2)**2/2.-phip(:,2),1,1)
      call color(4); call polyline(z,vz**2/2.-phip,imax)
!      call polymark(z,vz**2/2.-phip,1,1)
      call color(15)
      call pltend

    end subroutine doplots2
!***********************************************************************
    subroutine doplots(imax,wp0,w0,pt0)
      integer, dimension(2) :: imax
      if(lp)then; call pfset(3); else; call pfset(-3);
      endif
      
      call multiframe(3,1,0)
      call autoplot(z,r,imax)
      call axlabels('z','r')
!      call pltend
      
      call autoplot(z,vz**2/2.-phip,imax)
      call axlabels('z','v!dz!d!u2!u/2-!Af!@')
!      call pltend

      call autoplot(z,psi/cosh(z/4.)**4,imax)
      call pltend

      call multiframe(4,1,0)
      call autoplot(tv,vz**2/2.-phip,imax)
      call axlabels('t','v!dz!d!u2!u/2-!Af!@')
!      call pltend

      call autoplot(tv,vz**2/2.-psi/cosh(z/4.)**4,imax)
      call axlabels('t','v!dz!d!u2!u/2-!Af!@(r=r0)')
!      call pltend
      
      call autoplot(tv,(w-w0),imax)
      call axlabels('t','Dw')
!      call pltend
      
      call autoplot(tv,(pt-pt0)/pt0,imax)
      call axlabels('t','Dpt')
      call pltend

      call multiframe(0,0,0)
     
!      call autoplot(r,(vr**2+vt**2)/2.-phip,imax)
!      call axlabels('r','(v!dr!d!u2!u+v!d!At!@!d!u2!u)/2-!Af!@')
!      call pltend
     
!      call autoplot(z,Ep,imax)
!      call axlabels('z','Ep')
!      call pltend

      call charsize(.018,.018)
      call autoplot(xplot,yplot,imax)
      call axlabels('x','y')
      call pltend

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
end module orbitpoincare


!***********************************************************************
!****************************************************************************
   SUBROUTINE RKADVC(DX,X0,ND,YR1,icase,YR2,cond,IERR)
     use orbitpoincare
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
!***********************************************************************
subroutine orbitp
  use orbitpoincare
  character*78 string
  integer, dimension(2) :: imax

  wpt=0.
! We enter this point with w0,Bsqpsi,r,th,z set, vth=0 
  B=Bsqpsi*sqrt(psi)
  dt=0.047/B
!  dt=.01/B
  if(lprint)write(*,73)'w0=',w0,' Eropsi=',Eropsi,' Bsqpsi=',Bsqpsi,&
       ' psi=',psi,' r0=',y0(1) !,' rho0=', y0(4)/B

  if(idoplots.ne.0)then
  if(lp)then
     call pfset(3)    ! Display and output plots
  else
     call pfset(-3)   ! Run continuously outputing plots
  endif
  call pltinit(-.99999,1.,-1.,0.)
  call charsize(.018,.018)
! Thermal rho=1/B.  
     write(string,'(a,f4.2,a,f5.3,a,f4.2,a,f4.2,a,f4.2,a)') &
       'W=',w0,' E!dr!d/!Ay!@=',Eropsi &
       ,' !AW!@/!A)y!@=',Bsqpsi,' !Ar!@/r!d0!d=',1./(B*y0(1)) &
       ,' !Ay!@=',psi
  call boxtitle(string)
  call axis
  call axis2
  call axlabels('!Ac!@/!Ap!@','W!d!A|!@!d/!Ay!@')
endif

  do k=1,nwp
     t=0.
     if(nwp.le.2)then
        wp0=-(2*k-1)*wpf*psi
     else
        wp0=-k*psi/(nwp+1.)
        wp0=-psi*(k/(nwp+1.))**alpha
     endif
! Initial position (zero velocity)
     y1=y0
! Initial r/theta velocities
     vmod0=sqrt(2.*(w0-wp0))
     y1(4)=cos(vangle0*3.1415926/180.)*vmod0
     y1(5)=sin(vangle0*3.1415926/180.)*vmod0
     pt0=y1(1)*y1(5)-B*y1(1)**2/2.
! The following is assuming y0(5)=0. which is not necessarily true
! It needs to be the gyrocenter.
!  r0=y0(1)
! Instead we better approximate it as the place where y0(5) is zero
! although that is not perfect
     r0=sqrt(-2.*pt0/B)
     if(lprint.and.k.eq.nwp)write(*,*)'rc0=',r0
     phip1=psiofrz(r0,y0(3))
     y1(6)=sqrt(2.*(wp0+phip1))
! Obsolete:
! Set vr and vz consistent with Wp=-psi/2.
!     phip1=psiofrz(y0(1),y0(3))
!     y0(6)=sqrt(phip1)  ! Make wp=vz^2/2-phi equal to -phi/2
!     y0(4)=sqrt(2.*(phip1+w0)-y0(6)**2) ! Give vr the rest of energy.
! These velocities are changed in the energy scans below.  
!  w0=(y0(4)**2+y0(5)**2+y0(6)**2)/2.-phip1 ! The unnormalized W (not w)

     if(lprint)write(*,'(i3,a,f8.4,a,f8.4,$)')k,' Wp0=',wp0,' vz=',y1(6)
     j=0
     wpp=wp0
     do i=1,nstep
        if(i.le.nplot.and.k.le.2)then
           r(i,k) =y1(1)
           th(i,k)=y1(2)
           z(i,k) =y1(3)
           xplot(i,k)=r(i,k)*cos(th(i,k))
           yplot(i,k)=r(i,k)*sin(th(i,k))
           vr(i,k)=y1(4)
           vt(i,k)=y1(5)
           vz(i,k)=y1(6)
           phip(i,k)=psiofrz(y1(1),y1(3))
           Ep(i,k)=Eropsi*phip(i,k)
           w(i,k)=(vr(i,k)**2+vt(i,k)**2+vz(i,k)**2)/2.-phip(i,k)
           pt(i,k)=r(i,k)*vt(i,k)-B*r(i,k)**2/2.
           tv(i,k)=t
           imax(k)=i
        endif
        call RKADVC(dt,t,nvec,y1,icase,y2,cond,IERR)
        phipi=phip1
        wpi=y2(6)**2/2.-phip1 ! Parallel energy at r0
        if(wpi*wpp.lt.0.)then
           if(lprint)write(*,'(a,i6,a,f8.4,a,i4,a)') &
           ' Step=',i,' Wp sign flip' ,wpi,' after',j,' bounces'
           wpt=-wp0/psi
           wpp=wpi
!           exit
        endif
        if((.not.abs(y2(3)).lt.20.).or.(.not.abs(y2(1)).lt.2.*y0(1)) &
             .or..not.wpi.lt.10.)then
           if(lprint)write(*,'(a,i6,a,f6.2,a,f8.4,a,i4,a)')&
           ' Step=',i,' z=',y2(3), &
           ' Wp=',wpi,' after',j,' bounces'
           wpt=-wp0/psi
           exit
        endif
        if(y2(3)*y1(3).le.0)then ! We crossed the z=0 center.
           if(j.eq.nbouncemax)exit
           j=j+1
           phip1=psiofrz(r0,y2(3)) ! Evaluate new phi at ~gyrocenter
           f1=abs(y2(3))/(abs(y1(3))+abs(y2(3)))
           f2=abs(y1(3))/(abs(y1(3))+abs(y2(3)))
           tc(j)=t
           wp(j)= f1*(y1(6)**2/2.-phipi) &  ! old
                 +f2*(y2(6)**2/2.-phip1)    ! new
           xi(j)=atan2(y1(5),y1(4))
        endif
        t=t+dt
        y1=y2
     enddo
     if((i.eq.nstep+1.or.j.eq.nbouncemax).and.lprint)then
        write(*,'(a,i4,$)')' j=',j
        write(*,73)' tb/2=',t/(j-1.), &
             ' tcyc=',2.*3.1415926/B,' ratio=',t/(j-1.)*B/2/3.1415926
     endif
     if(idoplots.ne.0)call color(mod(k-1,14)+1)
     if(idoplots.ne.0)call polymark(xi/3.1415926,wp/psi,j-2,3)
  enddo
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
  if(idoplots.ne.0)call pltend
  
  if(nwp.le.2.and.idoplots.ge.2)call doplots2(imax,wp0,w0,pt0)  
  if(nwp.le.2.and.idoplots.ge.1)call doplots(imax,wp0,w0,pt0)  

!  583 format(i5,a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3)
  73 format(a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3)
end subroutine orbitp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine orbitmin
  ! Find the threshold confinement Omega for a psi range and Eropsi range.
  ! And plot it.
  use orbitpoincare
  integer, parameter :: nbscan=10,nbiter=5
  real :: OmegaLmax
  real, dimension(npmmax) :: xioL,opp,ppp
  character*40 thestr,thewpf
  character*2 cel

  cel='L '
  OmegaLmax=4.
  wpf=wpm
  call fwrite(wpf,iwidth,2,thewpf)
  nwp=1
  iwritetype=2         ! No writing from orbitp
  idoplots=0
  elmin=1./Eropsi
  do k=1,nerp
     Eropsi=1./(elmin*max(k-1.,0.5))
     do i=1,npm
        psi=psim*(i-.8)/(npm-.8)
        bmax=max(OmegaLmax*Eropsi,sqrt(psim))
2        do j=1,nbscan     ! Scan to find first confined orbit
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
        write(*,'(i3,a,f8.4,a,f8.4,a,f8.4)')i,' Eropsi=',Eropsi,' psi=',psi,' Omega=',bb
        psiplot(i,k)=psi
        omegaplot(i,k,2)=bb
        omegaplot(i,k,1)=bb/Eropsi
     enddo
  enddo
! Plotting
  call pfset(3)
  do i=1,2
     call pltinit(0.,psim,0.,OmegaLmax/i)
     call charsize(.018,.018)
     call axis
     call axis2
     call axlabels('!Ay!@', &
       cel(i:i)//'!AW!@ for W!d!A|!@t!d=-'//thewpf(1:iwidth)//'!Ay!@')
     call winset(.true.)
     do k=1,nerp
        Eropsi=1./(elmin*max(k-1.,0.5))
        call fwrite(1./Eropsi,iwidth,1,thestr)
        call color(k)
        call labeline(psiplot,omegaplot(:,k,i),npm,thestr,iwidth)
!        write(*,*)k,i,omegaplot(1:3,k,i)
     enddo
     call color(15)
     call dashset(1)
     opp=1.6/npm*[(j,j=1,npmmax)]
     if(i.eq.1)then
        xioL=(1+sqrt(1+2.*1/opp**2))
        ppp=opp**2*xioL*exp(-xioL)
        yl=0.1
        call legendline(0.5,yl-0.05,0,' Disruption')
     else
        ppp=opp**2
        yl=0.95
        call legendline(0.5,yl-0.05,0,' !AW!@=!Ay!@!u1/2!u')
     endif
     call legendline(0.5,yl,257,' Labels: L')
     call polyline(ppp,opp,npm)
     call dashset(0)
     call pltend
  enddo
  Eropsi=1./elmin
end subroutine orbitmin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program mainpoincare
  use orbitpoincare
  character*30 :: string
  real :: wn0=-10. ! Unset normalized w.
  do i=1,iargc()   ! Parse cmdline arguments.
     call getarg(i,string)
     if(string(1:2).eq.'-b')read(string(3:),*)Bsqpsi
     if(string(1:2).eq.'-O')then
        read(string(3:),*)Omega
        Bsqpsi=Omega/sqrt(psi)
        iwritetype=1
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
     if(string(1:2).eq.'-a')read(string(3:),*)alpha
     if(string(1:2).eq.'-v')read(string(3:),*)vangle0
     if(string(1:2).eq.'-c')lp=.not.lp
     if(string(1:2).eq.'-L')lo=.not.lo
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
     call orbitmin
  endif
  call exit
  
201 continue
  write(*,*)'Usage orbitpoincare [flags]'
  write(*,*)'-b... B/sqpsi, -p... psi, -q... sqrt(psi), -E... Er/psi, -r... r0'
  write(*,*)'-W... W0 (total energy), -w ... w0 (normalized alternate)'
  write(*,*)'-O... Omega, resets actual B/sqpsi value; use after all -p'
  write(*,*)'-nw... N-energies [if 1 plot orbits], -f... the single energy'
  write(*,*)'-a... set alpha power to concentrate wp values near zero. '
  write(*,*)'-v... set vangle0 degrees of initial velocity to direction r/y.'
  write(*,*)'-ns... N-steps, -c run continuously, -d do not write, -h help'

  write(*,*)'-mb[<wpm>,<psim>,<npm>] Find the minimum B that confines the orbit'
  write(*,*)'      of fractional depth wpm, up to psimax [0.5], in npm steps [100]'
  write(*,*)'-mE   set number of different Eropsi for -mb call.'
  write(*,*)'-L    toggle plotting disruption/resonance scaling'
end program mainpoincare

