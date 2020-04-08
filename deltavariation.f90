real function deltaofomega(Omega)
  real Omega  ! Omega/sqrt(psi) scaled
  deltaofomega=1/sqrt(Omega*exp(omega)+3.1415926/8.)
end function deltaofomega

real function powerform(Omega,n)
  real Omega  ! Omega/sqrt(psi) scaled
  integer n
  if(2.*Omega/n.lt.1.)then
     powerform=1/sqrt(Omega/(1.-2.*Omega/n)**(n/2)+3.1415926/8.)
  else
     powerform=1.e-10
  endif
end function powerform

program plotdeltaofomega
  integer, parameter :: np=500
  real, dimension(np) :: Omega,dofo,p2,p4,p6,p8,p10,P20

  do i=1,np
     Omega(i)=5.*(i-1.)/(np-1.)
     dofo(i)=deltaofomega(Omega(i))
     p2(i)=powerform(Omega(i),2)
     p4(i)=powerform(Omega(i),4)
     p6(i)=powerform(Omega(i),6)
     p8(i)=powerform(Omega(i),8)
     p10(i)=powerform(Omega(i),10)     
     p20(i)=powerform(Omega(i),20)     
  enddo

  !  call autoplot(Omega,dofo,np)
  call autoplot(Omega,dofo,np,.false.,.true.)
  call axlabels('!AW!@/!A)y!@','!Ad!@')
  call winset(.true.)
  call dashset(2)
  call polyline(Omega,p2,np,'2',1)
  call polyline(Omega,p4,np,'4',1)
  call polyline(Omega,p6,np,'6',1)
  call polyline(Omega,p8,np)
!  call polyline(Omega,p10,np)
  call polyline(Omega,p20,np)
  call pltend

end program plotdeltaofomega
