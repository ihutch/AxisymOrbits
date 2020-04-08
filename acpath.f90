!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Alternative module for using acpathdoc. Copy to code directory. Compile.
! Do 'use acpath' in f90 code. To turn on and off do nacpdo=0,1
! Cause acpath.o to be linked with your code before library.
module acpath
  integer, parameter :: nacpathmax=5000
  real, dimension(nacpathmax) :: xacpath,yacpath
  real :: cvacpath
  integer :: nacp,nacpdo=0   !default off
end module acpath
  ! replacement for accis dummy
subroutine acpathdoc(imax,xc,yc,cv)
  use acpath
  integer imax
  real xc(imax),yc(imax),cv
  if(nacpdo.eq.0)return
  if(imax.le.nacpathmax)then
     nacp=imax
     do i=1,imax
        xacpath(i)=xc(i)
        yacpath(i)=yc(i)
     enddo
     cvacpath=cv
  else
     write(*,*)'Not copied. Acpath length exceeded',imax,nacpathmax
  endif
end subroutine acpathdoc
