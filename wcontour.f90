  ! Plotting contours of Wparallel expressions. Old version1
module wcontour
  integer, parameter :: nxi=201,nwp=201,ncont=40
  real, dimension(nxi) :: xiv
  real, dimension(nwp) :: wpv
  real, dimension(nxi,nwp) :: ftotal,work
  real :: Er,B,w,wpr
  real, dimension(ncont) :: zclv
  
contains
  real function deepw(wpa,w,wpr)
    ! return the value of the deep trapping LHS with parameters
    ! wpa=-Wparallel/\psi; w=W/psi; wpr=Wparallel_resonant/psi 
    real wpa,w,wpr
    swpr=sqrt(wpr)
    deepw=atanh(sqrt(wpa*(w+1)/(w+wpa)))
    deepw=deepw-sqrt(w+1)*alog(sqrt(w+wpa)+sqrt(wpa))
    deepw=deepw-swpr*atanh(sqrt((w+wpa)/(w+1)))
    deepw=deepw*4./sqrt(w+1)
  end function deepw

  real function shallow(wpa,w,wpr) ! shallow trapping LHS
    real wpa,w,wpr
    swpr=sqrt(wpr)
    shallow=2.*(sqrt(w+wpa)-swpr*alog(sqrt(w+wpa)+sqrt(wpa)))
    shallow=shallow*2*3.1415926/8.
  end function shallow

  subroutine arraycalc
    do j=1,nwp
       wp=(j-1.)/(nwp-1.)*.995
       wpv(j)=wp
       dw=deepw(wp,w,wpr)+shallow(wp,w,wpr)
!       dw=+shallow(wp,w,wpr)
!       dw=+deepw(wp,w,wpr)
       do i=1,nxi
          xiv(i)=(-1.+2.*(i-1.)/(nxi-1.))
          xi=3.1415926*xiv(i)
          ftotal(i,j)=dw-Er*sin(xi)
       enddo
    enddo
  end subroutine arraycalc  
end module wcontour

program docontours
  use wcontour

  wpr=.8
  w=3
  Er=.01
 
  call arraycalc
  if(nxi.le.10)then
     write(*,101)xiv
     write(*,101)wpv
     write(*,101)ftotal
  endif
  call pfset(3)
  call pltinit(-1.,1.,-1.,0.)

  do k=1,3
     if(k.eq.1)wpr=.9
     if(k.eq.2)wpr=.5
     if(k.eq.3)wpr=.05
     call arraycalc
     icsw=1
     icl=11
     call color(k+1)
     do i=1,icl
        zclv(i)=ftotal(1,nint(wpr*nwp))+ Er*(-1.+4.*float(i-icl)/(1-icl))
     enddo
     icl=-icl
     call contourl(ftotal,work,nxi,nxi,nwp,zclv,icl,xiv,-wpv,icsw) 
  enddo

  call color(15)
  call axis
  call axis2
  call axlabels('!Ac!@/2!Ap!@','W!d!A|!@!d/!Ay!@ (= -w!d!a|!@!d)')
  call pltend
  
101 format(10f8.4)
end program docontours
