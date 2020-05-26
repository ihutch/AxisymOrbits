  ! Code to plot the domain in b and perturbation strength the W_\parallel
  ! boundary value between trapped and detrapped orbits.

module trapdomain
  integer, parameter :: nb=300, nE=300
  real, dimension(nb) :: b                ! Normalized to psi^(1/2)
  real, dimension(nE) :: Erv              ! Normalized to psi^(3/2)
  real, dimension(nE,nb) :: wpbdy,work
  integer, dimension(nE,nb) :: nlowres
  real, dimension(nb) :: wpbb,wpbr
  integer, dimension(nb) :: nbb
  real, dimension(nE) :: vperp(nE),vpar(nE)
  real :: bmax=3,Ervmax=.4,vemax=2.4,psi=1.,Eropsi=.1
  integer :: nres,ipfset=3,npsi=5
  character*30 :: filename
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine findwpar(Ervv,bv,wpbdv,nresv,swprv)
    ! Given Ervv,bv, find wpbdv (wp boundary),
    ! nresv (lowest overlapped resonance No), swprv (sqrt(resonance energy)) 
    swpr=2.  ! prevent exit on first cycle 
    deltan=0.
    ff=1.
    do n=2,40,2
       swprp=swpr
       deltanp=deltan
       swpr=sqrt(((2.*bv**2/n**2)**(-0.25)+1-2.**.25)**(-4))
       deltan=sqrt(Ervv*2./n)/sqrt(swpr*(n/2.)**ff/(1-swpr)**(n/2)+3.1415926/8.)
       if(swpr.lt.1)then
          if(swprp-deltanp.lt.swpr+deltan)then
             exit
          endif
       endif
    enddo
    nresv=n-2
    wpbdv=-min(1.0,(swprp+deltanp)**2)
    swprv=swprp
  end subroutine findwpar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine doWcontours
    integer, parameter :: ncl=20
    real, dimension(ncl) :: zclv
    do j=1,nb
       b(j)=bmax*float(j)/nb
       do i=1,nE
          Erv(i)=Ervmax*float(i)/nE
          call findwpar(Erv(i),b(j),wpbdy(i,j),nlowres(i,j),swprv)
       enddo
    enddo
    if(nE.le.10)then
       write(*,14)nlowres
       write(*,84)wpbdy
    endif
    
!    call accisinit ! Not needed when pfset has been called.
    call charsize(.018,.018)
    call accisgradinit(-32000,0,-65000,97000,150000,65500)
    call autocolcont(wpbdy,nE,nE,nb)
    call legendline(-.34,.7,258,'W!d!A|!Bt!@!d/!Ay!@')
    icl=0
    zclv(1)=5.
    icsw=1
    call scalewn(0.,Ervmax,0.,bmax,.false.,.false.)
    call axlabels('E!dr0!dv!d!A`!@!d/!Ay!@!u3/2!u=(v!d!A`!@!d/!A)y!@)/L!d!A`!@!d','')
    call legendline(-.15,.58,258,'!AW!@/!A)y!@')
    call axis
    call axis2
    call color(9)
    call contourl(wpbdy,work,nE,nE,nb,zclv,icl,Erv,b,icsw)
    call color(13)
    call vecw(0.,0.,0)
    call vecw(Ervmax,Ervmax,1)
    call legendline(.68,.05,258,'Invalid: !Ar!@>L!d!A`!@!d')
    call pltend()
    
    
14   format(10i4)
84  format(10f8.4)
  end subroutine doWcontours
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine lineout(nfile)
    integer, parameter :: nchmax=30
    real, dimension(nchmax) :: bcheck,wpt
    logical :: linit=.false.
    character*10 :: string
    
    j=0
    open(22,file=filename,status='old',err=101)
    do j=1,nchmax
       read(22,*,end=100,err=100)bcheck(j),wpt(j),Erovpsi3
    enddo
100 j=j-1
    Ervmax=Erovpsi3
    write(*,'(3a,i3,a,f7.3)')'Read ',filename,' Lines:',j,'  Ervmax=',Ervmax
101 continue
    bmin=Ervmax
    do i=1,nb
       b(i)=bmin+(bmax-bmin)*float(i)/nb
       call findwpar(Ervmax,b(i),wpbb(i),nbb(i),swprv)
       wpbr(i)=-swprv**2
    enddo
    if(.not.linit)then
       call charsize(.018,.018)
       call pltinit(0.,b(nb),-1.,0.)
       call axis
       call axis2
       call axlabels('!AW!@/!A)y!@','W!d!A|!Bt!@!d/!Ay!@')
       linit=.true.
       call legendline(0.6,0.36,258, ' E!dr0!dv!d!A`!@!d/!Ay!@!u3/2!u=')
    endif
    call color(nfile)
    call polyline(b,wpbb,nb)
    if(j.gt.0)then
       call dashset(4)
       call polyline(bcheck,-wpt,j)
       write(string,'(f4.2)')Ervmax
       call legendline(0.7,0.36-nfile*0.08,-nfile,' '//string)
       call dashset(0)
       call polymark(bcheck,-wpt,j,nfile)
!       call polyline(b,wpbr,nb)
    endif
    
  end subroutine lineout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine boundaryv(nv)
    character*100 string
    ! Plot the trapped passing boundary on vpar/vperp domain.
    if(nv.gt.nE)stop 'nv>nE error'
    call pltinit(0.,1.,0.,vemax)
    call charsize(0.018,0.018)
    call axis
    call axlabels( &
         'v!d!A|!@0!Bt!@!d/!A)!@(2!Ay!@) = (1+W!d!A|!Bt!@!d/!Ay!@)!u1/2!u' &
         ,'')   ! = !ArW!@
    call legendline(-.1,.55,258,'v!d!A`!@!d')
    ! The vperp at which rho=L so vperp=B/Eropsi
    verhoL=bmax/Eropsi
    rhomaxL=vemax/bmax*Eropsi
    write(string,'(a,f4.1,a,f4.2,a,f4.2,a)')'L!d!A`!@!d=',1./Eropsi &
         ,' !AW!@=',bmax,'   Labels: !Ay!@'
    call axptset(1.,0.)
    call ticrev
    call altyaxis(1/bmax,1/bmax)
    call legendline(1.03,.55,258,'!Ar!@')
    call ticrev
    call axptset(0.,0.)
    call termchar(string)
    call boxtitle(string)
    psimax=psi
    do j=1,npsi
       psi=psimax*j/npsi
       write(*,'(a,4f8.4)')'psi,bmax,Eropsi,rhomaxL',psi,bmax,Eropsi,rhomaxL
       do i=1,nv
          ve=vemax*i/float(nv)
          Erv(i)=Eropsi*ve/sqrt(psi) 
          bv=bmax/sqrt(psi)
          call findwpar(Erv(i),bv,wpbdv,nresv,swprv)
          vperp(nv+1-i)=ve
          vpar(nv+1-i)=sqrt(wpbdv+1)
          if(nv.le.10)write(*,'(i3,3f8.4,i3)')i,vpar(i),vperp(i),wpbdv,nresv
       enddo
       call color(j)
       call dashset(j)
       write(string,'(f4.2)')psi
       call termchar(string)
       call labeline(vpar,vperp,nv,string,4)
!       call polyline(vpar,vperp,nv)
    enddo
    call pltend
  end subroutine boundaryv

  
end module trapdomain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program trapmain
  use trapdomain
  character*30 string
  integer :: ipno=1
  filename=' '
  nfile=0
  do i=1,iargc()   ! Parse cmdline arguments.
     call getarg(i,string)
     if(string(1:2).eq.'-b')read(string(3:),*)bmax
     if(string(1:2).eq.'-E')read(string(3:),*)Ervmax
     if(string(1:2).eq.'-L')read(string(3:),*)perpL
     if(string(1:2).eq.'-c')ipfset=-3
     if(string(1:2).eq.'-t')read(string(3:),*)ipno
     if(string(1:2).eq.'-h')goto 201
     if(string(1:2).eq.'-?')goto 201
     if(string(1:1).ne.'-')then
        call pfset(ipfset) ! Have to set here before first file call.
        read(string,*,err=101)filename
        ipno=2
        nfile=nfile+1
        call lineout(nfile)
101     continue
     endif
  enddo
  Eropsi=1./perpL
  call pfset(ipfset)

  if(ipno.eq.2)then
     if(filename(1:1).eq.' ')call lineout(15)
     call pltend
  elseif(ipno.eq.1)then
     call doWcontours
  elseif(ipno.eq.3)then
     nv=100
     call boundaryv(nv)
  endif


  
  call exit
201 write(*,*)'-b... bmax, -E... Ervmax -t... plot type'
end program trapmain
