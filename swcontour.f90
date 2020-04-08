! Plotting contours of Wparallel expressions.
! Based on better interpolation but more complicated functions.
! And giving the harmonics, by iteration of a recursion relation.
module swcontour
  integer, parameter :: nxi=201,nwp=201,ncont=40
  real, dimension(nxi) :: xiv
  real, dimension(nwp) :: wpv,Fnarray,Fplot
  real, dimension(nxi,nwp) :: ftotal,work
  real :: Er,B,w,wpr
  real, dimension(ncont) :: zclv
  integer :: n ! Harmonic number 2,4,6 ...
  
contains
  real function h1(w,wp)
    h1=alog(sqrt(w+wp)+sqrt(wp))
  end function h1
  real function h2(w,wp)
    h2=alog(1-sqrt(wp))-alog(sqrt(w+1)*sqrt(w+wp)+w+sqrt(wp))
  end function h2
  
  real function Fn(n,w,wp,wpr)
! Obsolete form based on Wolfram results. Good for testing.
    if(n.eq.2)then      !g2
       gn=-2.*sqrt(w+wp)+2.*(sqrt(wpr)-1.)&
            *(h1(w,wp)+h2(w,wp)/sqrt(w+1))
    elseif(n.eq.4)then  !g4
       gn=-2.*(sqrt(wpr)-1.)*sqrt(w+wp)/((w+1)*(1.-sqrt(wp)))  &
            -(2.*w*(sqrt(wpr)-2.)-2.)*h2(w,wp)/(w+1)**1.5    &
            +2.*h1(w,wp)
    elseif(n.eq.6)then  !g6
       gn=( w*(  sqrt(wpr)*(1-2.*sqrt(wp))+4.*sqrt(wp)-3.  )  &
            +sqrt(wpr)*(sqrt(wp)-2)+sqrt(wp) )*sqrt(w+wp)       &
            /((w+1)*(1.-sqrt(wp)))**2 &
            -w*(2.*w+3.*sqrt(wpr)-1)*h2(w,wp)/(w+1)**2.5
    elseif(n.eq.8)then
       gn=(-(1-sqrt(wp))**2*( -6.*w**2+w*(10.-13.*sqrt(wpr))+2.*sqrt(wpr)+1.) &
            +(w+1.)*(1.-sqrt(wp))*(3.*w*(sqrt(wpr)-2.)-2.*sqrt(wpr)-1) &
            +2.*(w+1.)**2*(1-sqrt(wpr)) )*sqrt(w+wp) &
            /(3.*(w+1)**3*(1-sqrt(wp))**3)           &
            +w*(w*(sqrt(wpr)-4.)-4.*sqrt(wpr)+1)*h2(w,wp)/(w+1)**3.5
    else
       write(*,*)'Unsupported harmonic number',n
       stop
    endif
    Fn=n*3.1415926/(16.*sqrt(2.))*2*(sqrt(w+wp)-sqrt(wpr)*h1(w,wp)) !g0
!        Fn=Fn+n**2/(4.*sqrt(2.))*gn ! Old error.
    Fn=Fn+n**2/(2.*sqrt(2.))*gn
  end function Fn

  real function Fniter(n,w,wp,wpr)
    ! version to calculate by iteration.  Working and agrees.
    real, dimension(0:20) :: gi,gj
    if(n.lt.2.or.n.gt.40)then
       write(*,*)'Harmonic choice n=',n,'out of range'
       stop
    endif
    gj0=sqrt(w+wp) 
    gi(0)=alog(gj0+sqrt(wp))  
    gi(1)=-(alog(1-sqrt(wp))-alog(sqrt(w+1)*sqrt(w+wp)+w+sqrt(wp)))/sqrt(w+1)
    gj(0)=gj0
    gj(1)=gi(1)-gi(0)
    do m=1,n/2-1
       gi(m+1)=( gj0/(1-sqrt(wp))**m+(2.*m-1)*gi(m)-(m-1)*gi(m-1))/(m*(w+1) )
       gj(m+1)=gi(m+1)-gi(m)
    enddo
    gn=2.*( (1-sqrt(wpr))*gj(m)-gj(m-1) )

    Fniter=n*3.1415926/(16.*sqrt(2.))*2*(sqrt(w+wp)-sqrt(wpr)*h1(w,wp)) !g0
    Fniter=Fniter+n**2/(2.*sqrt(2.))*gn
  end function Fniter

  
  subroutine arraycalc
    do j=1,nwp
       wp=(j-1.)/(nwp-1.)*.995
       wpv(j)=wp
       dw=Fniter(n,w,wp,wpr)
       Fnarray(j)=dw
       if(.false.)then  ! Fn comparison only
          dwi=Fn(n,w,wp,wpr)
          write(*,*)j,dw,dwi,dw-dwi
          if(j.eq.nwp)stop
       endif
       do i=1,nxi
          xiv(i)=(-1.+2.*(i-1.)/(nxi-1.))
          xi=3.1415926*xiv(i)
          ftotal(i,j)=dw-Er*sin(xi)
       enddo
    enddo
  end subroutine arraycalc
  
end module swcontour

program docontours
  use swcontour
  use acpath
  character*100 string
  real, dimension(4) :: wprs
  logical :: lp=.true.
  data wprs/.2,.5,.05,.9/

  ! Defaults
  wp2=0.
  nFs=0
  nn=6
  nwr=0
  bdef=.8
  w=1
  Er=.01
  b=sqrt(bdef)
  do i=1,iargc()  ! Parse cmdline arguments.
     call getarg(i,string)
     if(string(1:2).eq.'-b')read(string(3:),*)b
     if(string(1:2).eq.'-w')read(string(3:),*)w
     if(string(1:2).eq.'-E')read(string(3:),*)Er
     if(string(1:2).eq.'-R')read(string(3:),*)wp2
     if(string(1:2).eq.'-n')read(string(3:),*)nn
     if(string(1:2).eq.'-m')read(string(3:),*)nwr
     if(string(1:2).eq.'-F')read(string(3:),*)nFs
     if(string(1:2).eq.'-c')lp=.not.lp
     if(string(1:2).eq.'-h')goto 201
  enddo
  if(lp)then
     call pfset(3)
  else
     call pfset(-3)
  endif
  wp2=b**2 ! Resonance n is n*omegab=Omega (fixed): wpn=(2/n)^2*wp2
  write(*,84)'b=',B,' wp2=',wp2,' Er=',Er
  84 format(a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4)
  if(nwr.ne.0)then
     call pltinit(-.99999,1.,-1.,0.)
     call charsize(0.018,0.018)
     n=2
     do k=1,min(nwr,4)  ! n=2 cases with different wpr
        wpr=wprs(k)
        crange=4
        if(k.eq.1)call jdrwstr(wx2nx(1.03),wy2ny(0.02),'w!d!A|!@R!d=',1.)
        call fwrite(wpr,iwidth,2,string)
        call color(k)
        call arraycalc
        icsw=1
        icl=11
        nf=nint(wpr*nwp)+1
        do i=1,icl
           zclv(i)=ftotal(1,nf) &
                + Er*(-1.+crange*float(i-icl)/(1-icl))
        enddo
        icl=-icl
        call contourl(ftotal,work,nxi,nxi,nwp,zclv,icl,xiv,-wpv,icsw)
        call jdrwstr(wx2nx(1.07),wy2ny(-wpr),string(1:lentrim(string)),.5)
        if(k.eq.nFs)then
           Fplot=Fnarray-minval(Fnarray)
           write(*,'(a,i3,a,f8.4,a,f8.4)') &
                'Plotting Fn for n=',n,' wpr=',wpr,' w=',w
        endif
     enddo
     call color(15)
     call axis
     call axis2
     call axlabels('!Ac!@/!Ap!@','W!d!A|!@!d/!Ay!@')
     call pltend
  endif

  call pltinit(-.99999,1.,-1.,0.)
  call charsize(0.018,0.018)
  write(string,'(a,f5.2,a,f5.3,a,f5.3,a,f5.3)') &
  'W/!Ay!@=',w,' E!dr!d/!Ay!@=',Er &
       ,' !AW!@/!A)y!@=',B
  call boxtitle(string(1:lentrim(string)))

  wprm=0.
  do k=1,nn ! Fixed b cases.
     n=2*k
!     wpr=wp2*4/n**2
     wpr=((2.*B**2/n**2)**(-0.25)+1-2.**.25)**(-4)
485 format(a,i4,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f8.5)
     write(*,485)'n=',n,' wpr=',wpr,' 4B^2/n^2=',4.*B**2/n**2
     if(abs(wpr-wprm).lt.0.03)exit
     wprm=wpr
     crange=4
     if(k.eq.1)call jdrwstr(wx2nx(1.07),wy2ny(0.02),'n=',0.)
     call iwrite(n,iwidth,string)
     nf=nint(wpr*nwp)+1
!     if(nf.ge.1 .and. nf.le.nwp)then  ! This resonance is in-range
     if(nf.ge.1.and.nf.lt.1.3*nwp)then
        if(nf.gt.nwp)nf=nwp
        call color(k)
        call arraycalc
        icsw=1
        icl=11
        do i=1,icl
           zclv(i)=ftotal(1,nf) &
                + Er*(-1.+crange*float(i-icl)/(1-icl))
        enddo
        icl=-icl
        call contourl(ftotal,work,nxi,nxi,nwp,zclv,icl,xiv,-wpv,icsw)
        call jdrwstr(wx2nx(1.07),wy2ny(-wpr),string(1:lentrim(string)),0.)
        icl=-1
        zclv(1)=ftotal(1,nf)+Er*0.98
        nacpdo=1  ! Turn on path following.
        call contourl(ftotal,work,nxi,nxi,nwp,zclv,icl,xiv,-wpv,icsw)
        nacpdo=0  ! Turn off
        write(*,'(a,i3,3f8.4)')'Sepx range, width',k,minval(yacpath(1:nacp)) &
             ,maxval(yacpath(1:nacp)) &
             ,maxval(yacpath(1:nacp))-minval(yacpath(1:nacp))
        zclv(2)=ftotal(1,nf)+Er*1.02
        call contourl(ftotal,work,nxi,nxi,nwp,zclv(2),icl,xiv,-wpv,icsw)
     endif
     if(k.eq.nFs)then
        Fplot=Fnarray-minval(Fnarray)
        write(*,'(a,i3,a,f8.4,a,f8.4)') &
             'Plotting Fn for n=',n,' wpr=',wpr,' w=',w
     endif
  enddo
  call color(15)
  call axis
  call axis2
  call axlabels('!Ac!@/!Ap!@','W!d!A|!@!d/!Ay!@')
  if(nn.gt.0)call pltend

  if(nFs.ne.0)then
     call pltinit(0.,1.,0.,1.)
     call axis
     !     call autoplot(wpv,Fplot,nwp)
     call axlabels('wp','Fn')
     call winset(.true.)
     call polyline(wpv,Fplot,nwp)
    call pltend()
  endif
  call exit
201 continue
  write(*,*)'Usage ./swcontour [-b..] [-w..] [-E..] [-R..] [-h] [...]'
  write(*,*)'Flags set -b B-field, -w total energy, Er/psi, -R resonant energy.'
  write(*,*)'-n #harmonics [default 6], -m #wprs, -F Function# to plot'
  write(*,*)'-c run/plot without stopping'
  !101 format(10f8.4)
end program docontours
