  ! Plot an analytic island in W vs \xi.

program simpleisland

  integer, parameter :: nxi=201,nc=10
  real, dimension(nxi) :: xi,Wp,Wm
  real, dimension(nc) :: Cn
  integer :: i,j
  character*10 string
  
  W0=-.5
  psi=1.
  EvWo=.02
  call pfset(3)
  call fwrite(EvWo,iwidth,3,string)
  call multiframe(1,2,1)
  xi=(/(-1.+2.*(j-1.)/(nxi-1.),j=1,nxi)/)  !*3.1415926
  call pltinit(-.99999,1.,-psi,0.)
  call charsize(.018,.018)
  call axis
  call axis2
  call axlabels('!Ac!@/!Ap!@','W!d!A|!@!d')
  call legendline(1.02,.5,258,'W!dR!d')
  call boxtitle(string)
  do i=1,nc
     Cn(i)=(-1.+2.*(i-1.)/(nc-5.))
     Wp=W0+sqrt(max(EvWo*(-sin(3.1415926*xi)+Cn(i)),0.))
     Wm=W0-sqrt(max(EvWo*(-sin(3.1415926*xi)+Cn(i)),0.))
!     Wm=2.*W0-Wp
     call polyline(xi,Wp,nxi)
     call polyline(xi,Wm,nxi)
  enddo
  write(*,'(10f8.3)')Cn
  call pltend
  
end program simpleisland
