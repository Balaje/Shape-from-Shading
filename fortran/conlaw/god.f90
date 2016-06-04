  ! Program to solve the Burger's equation
  ! u_t + f(u)_x = 0
  ! u(x,0) = { 0, x < 0
  !            1, x > 0
  !          }

program god

  use fun, only : f, fluxgod
  implicit none

  integer , parameter :: N = 100
  real(kind = 8), parameter :: a = -1., b = 2., eps = 1.
  real(kind = 8) :: h, delt, tf, t0, ul = 0., ur = 1.
  real(kind = 8), dimension(N-1) :: x, u, unew
  integer :: i,j,M,k


  h = (b-a)/N

  do i = 1,N-1
     x(i) = a + (i-1)*h
     if(x(i) < 0.) then
        u(i) = ul
     elseif(x(i) >= 0.) then
        u(i) = ur
     endif
  enddo

  delt = eps*h

  tf = 1.
  t0 = 0.

  M = (tf-t0)/delt

  do k = 1,M
     unew(1) = u(1) - eps*(fluxgod(f,u(1),u(2)) - fluxgod(f,ul,u(1)))
     do j=2,N-2
        unew(j) = u(j) - eps*(fluxgod(f,u(j),u(j+1)) - fluxgod(f,u(j-1),u(j)))
     enddo
     unew(N-1) = u(N-1) - eps*(fluxgod(f,u(N-1),ur) - fluxgod(f,u(N-2),u(N-1)))

     u = unew
  enddo

  open(unit = 1, file="sol_godunov.txt", status="unknown", action="write")
  do i = 1,N-1
     write(unit = 1, fmt=*)x(i), unew(i)
  enddo
end program god
