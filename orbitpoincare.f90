  ! Integrate orbits in a specified potential to generate a poicare plot
module orbitpoincare
  real :: Bsqpsi=.85, Eropsi=.01, psi=1. 
  real :: r0,wpf=0.5,B,w0=1.,alpha=1.
  logical :: lp=.true.,lprint=.true.
  integer, parameter :: nbouncemax=200  ! Max bounces to store.
  integer, parameter :: nplot=25000     ! Max steps to store.
  integer, parameter :: nvec=6
  integer :: nwp=39,nstep=1000000       ! Default maximum steps to take.
  real, dimension(nplot) :: r,th,z,vr,vt,vz,w,tv,pt,xplot,yplot,phip,Ep
  real, dimension(nbouncemax) :: wp,xi,tc
  real, dimension(nvec) :: y0,y1,y2
  data y0/10.,0.,0.,2.5,0.,0.9/ !Initial orbit conditions defaults
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
!      psir=psi*(1.+(r0-y(1))*Eropsi) ! Inconsistent with Er/psi=const.
      psir=psi*exp((r0-y(1))*Eropsi)  ! Local psi(r) 
      sc=1./cosh(y(3)/4.)
      phi=psir*sc**4
      E(1)=Eropsi*phi
      E(2)=0.
      E(3)=psir*sinh(y(3)/4.)*sc**5
    end subroutine getfield
!***********************************************************************
    real function psiofrz(rh,zh)
      cc=1./cosh(zh/4.)
      psir=psi*exp((r0-rh)*Eropsi)
      psiofrz=cc**4*psir
    end function psiofrz
!***********************************************************************
!***********************************************************************
    subroutine doplots(imax,wp0,w0,pt0)
      if(lp)then; call pfset(3); else; call pfset(-3); endif
      call autoplot(z,r,imax)
      call axlabels('z','r')
      call pltend
      
      call autoplot(z,vz**2/2.-phip,imax)
      call axlabels('z','v!dz!d!u2!u/2-!Af!@')
      call pltend
      call autoplot(tv,vz**2/2.-phip,imax)
      call axlabels('t','v!dz!d!u2!u/2-!Af!@')
      call pltend
     
      call autoplot(r,(vr**2+vt**2)/2.-phip,imax)
      call axlabels('r','(v!dr!d!u2!u+v!d!At!@!d!u2!u)/2-!Af!@')
      call pltend

      call autoplot(z,vz**2/2.-psi/cosh(z/4.)**4,imax)
      call axlabels('z','v!dz!d!u2!u/2-!Af!@(r=r0)')
      call pltend
      call autoplot(tv,vz**2/2.-psi/cosh(z/4.)**4,imax)
      call axlabels('t','v!dz!d!u2!u/2-!Af!@(r=r0)')
      call pltend
      
      call autoplot(tv,(w-w0),imax)
      call axlabels('t','Dw')
      call pltend

      call autoplot(tv,(pt-pt0)/pt0,imax)
      call axlabels('t','Dpt')
      call pltend

     
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
      call scale3(xmin,xmax,ymin,ymax,-5.,5.)
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
  real :: wpt=0.
  
! We enter this point with w0,Bsqpsi,r,th,z set, vth=0 
  B=Bsqpsi*sqrt(psi)
  dt=0.047/B
  pt0=y0(1)*y0(5)-B*y0(1)**2/2.
  r0=y0(1)
! Set vr and vz consistent with Wp=-psi/2.
  phip1=psiofrz(y0(1),y0(3))
  y0(6)=sqrt(phip1)  ! Make wp=vz^2/2-phi equal to -phi/2
  y0(4)=sqrt(2.*(phip1+w0)-y0(6)**2) ! Give vr the rest of energy.
! These velocities are changed in the energy scans below.  
!  w0=(y0(4)**2+y0(5)**2+y0(6)**2)/2.-phip1 ! The unnormalized W (not w)
 
  if(lprint)write(*,73)'w0=',w0,' Eropsi=',Eropsi,' Bsqpsi=',Bsqpsi,&
       ' psi=',psi,' r0=',r0, ' rho0=', y0(4)/B

! Thermal rho=1/B.  
  write(string,'(a,f4.2,a,f5.3,a,f4.2,a,f4.2,a,f4.2,a)') &
       'W=',w0,' E!dr!d/!Ay!@=',Eropsi &
       ,' !AW!@/!A)y!@=',Bsqpsi,' !Ar!@/r!d0!d=',1./(B*r0) &
       ,' !Ay!@=',psi
  if(lp)then
     call pfset(3)
  else
     call pfset(-3)
  endif
  call pltinit(-.99999,1.,-1.,0.)
  call charsize(.018,.018)
  call boxtitle(string)
  call axis
  call axis2
  call axlabels('!Ac!@/!Ap!@','W!d!A|!@!d/!Ay!@')

  
  do k=1,nwp
     t=0.
     y1=y0
     if(nwp.le.1)then
        wp0=-wpf*psi
     else
        wp0=-k*psi/(nwp+1.)
        wp0=-psi*(k/(nwp+1.))**alpha
     endif
     y1(6)=sqrt(2.*(wp0+phip1))
     y1(4)=sqrt(2.*(w0-wp0))
     if(lprint)write(*,'(i3,a,f8.4,a,f8.4,$)')k,' Wp0=',wp0,' vz=',y1(6)
     j=0
     do i=1,nstep
        if(i.le.nplot)then
           r(i) =y1(1)
           th(i)=y1(2)
           z(i) =y1(3)
           xplot(i)=r(i)*cos(th(i))
           yplot(i)=r(i)*sin(th(i))
           vr(i)=y1(4)
           vt(i)=y1(5)
           vz(i)=y1(6)
           phip(i)=psiofrz(y1(1),y1(3))
           Ep(i)=Eropsi*phip(i)
           w(i)=(vr(i)**2+vt(i)**2+vz(i)**2)/2.-phip(i)
           pt(i)=r(i)*vt(i)-B*r(i)**2/2.
           tv(i)=t
           imax=i
        endif
        phipi=phip1
        wpi=y1(6)**2/2.-phip1 ! Parallel energy at r0
        if(wpi.gt.0)then
           if(lprint)write(*,'(a,i6,a,f8.4,a,i5,a)')' Step=',i,' Positive Wp' &
                ,wpi,' after',j,' bounces'
           wpt=-wp0/psi
           exit
        endif
        if(abs(y1(3)).gt.20.)then
           if(lprint)write(*,'(a,i6,a,f6.2,a,f8.4,a,i4,a)')&
           ' Step=',i,' z=',y1(3), &
           ' Wp=',wpi,' after',j,' bounces'
           wpt=-wp0/psi
           exit
        endif
        call RKADVC(dt,t,nvec,y1,icase,y2,cond,IERR)
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
     call color(mod(k-1,14)+1)
     call polymark(xi/3.1415926,wp/psi,j-2,3)
  enddo
  write(*,'(3f7.4,a)')Bsqpsi,wpt,Eropsi*sqrt(2.*(w0/psi-wpt)) &
       ,'   b,wpt,Erovpsi3'


  call pltend
  
  if(nwp.eq.1)call doplots(imax,wp0,w0,pt0)  

!  583 format(i5,a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3)
  73 format(a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3)
end subroutine orbitp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program mainpoincare
  use orbitpoincare
  character*30 :: string
  real :: wn0=-10. ! Unset normalized w.
  do i=1,iargc()   ! Parse cmdline arguments.
     call getarg(i,string)
     if(string(1:2).eq.'-b')read(string(3:),*)Bsqpsi
     if(string(1:2).eq.'-p')read(string(3:),*)psi
     if(string(1:2).eq.'-E')read(string(3:),*)Eropsi
     if(string(1:2).eq.'-r')read(string(3:),*)y0(1)
     if(string(1:2).eq.'-f')read(string(3:),*)wpf
     if(string(1:2).eq.'-W')read(string(3:),*)W0
     if(string(1:2).eq.'-w')read(string(3:),*)wn0
     if(string(1:2).eq.'-a')read(string(3:),*)alpha
     if(string(1:3).eq.'-nw')read(string(4:),*)nwp
     if(string(1:3).eq.'-ns')read(string(4:),*)nstep
     if(string(1:3).eq.'-c')lp=.not.lp
     if(string(1:3).eq.'-d')lprint=.not.lprint
     if(string(1:2).eq.'-h')goto 201
     if(string(1:2).eq.'-?')goto 201
  enddo
  if(wn0.gt.-10.)w0=psi*wn0

  call orbitp
  call exit
  
201 continue
  write(*,*)'Usage orbitpoincare [flags]'
  write(*,*)'-b... B/sqpsi, -p... psi, -E... Er/psi, -r... r0, -W... W0, -w ... w0'
  write(*,*)'-nw... N-energies [if 1 plot orbits], -f... the single energy'
  write(*,*)'-ns... N-steps, -c run continuously -d do not write -h help'
  write(*,*)'-a... set alpha power to concentrate wp values near zero. '
end program mainpoincare

