  ! Program to solve the Burger's equation using Lax Friedrich scheme
  !    u_t + f(u)_x = 0
  !   u(x,0) = { 0, x<0
  !              1, x>0 }
  !

program lax

  use fun , only : f, fluxlax
  implicit none

  integer , parameter :: N = 100
  real (kind = 8), parameter :: a = -1., b = 2., eps = 1.
  real (kind = 8) :: delt, h, tf, t0, ul = 0., ur = 1.
  real (kind = 8), dimension(N-1) :: x,u,unew
  integer :: i,j,k,M

  h = (b-a)/N

  do i=1,N-1
     x(i) = a+(i-1)*h

     if(x(i)<0) then
        u(i) = ul
     elseif(x(i) >= 0) then
        u(i) = ur
     endif
  enddo

  tf = 1.
  t0 = 0.

  delt = eps*h
  M = (tf-t0)/delt

  do k = 1,M
     unew(1) = u(1) - eps*(fluxlax(f,u(1),u(2),eps) - fluxlax(f,ul,u(1),eps))
     do j = 2,N-2
        unew(j) = u(j) - eps*(fluxlax(f,u(j),u(j+1),eps) - fluxlax(f,u(j-1),u(j),eps))
     enddo
     unew(N-1) = u(N-1) - eps*(fluxlax(f,u(N-1),ur,eps) - fluxlax(f,u(N-2),u(N-1),eps))

     u = unew
  enddo

  open(unit = 2, file = "solution_lax.txt", status="unknown", action = "write")
  do i=1,N-1
     write(unit = 2, fmt = *)x(i), unew(i)
  enddo


end program lax
