! Program to solve the Burger's equation using Godunov's Scheme
!     u_t + f(u)_x = 0
!       u(x,0) = { 1, x < 0
!                  0, x >= 0 }
!

program godunov

  use functions, only : f
  use flux, only : fgod

  integer, parameter :: N = 500
  real(kind = 8) :: a = 0., b = 1., h, delt, tf, t0, r, lt = 1., rt = 0.
  real(kind = 8), dimension(N-1) :: v,vnew,x
  integer :: M, i, j

  open(unit = 1, file="god.txt", status="unknown", action="write")

  h = (b-a)/N

  ! CFL Condition
  delt = 1*h

  ! Defining the mesh and the initial condition
  pi = acos(-1.)
  do i = 1,N-1
    x(i) = a + i*h
    v(i) = sin(pi*x(i))
  enddo

  tf = 0.2
  t0 = 0.

  M = (tf-t0)/delt

  r = delt/h

  do j = 1,M
    vnew(1) = v(1) - r*(fgod(f,v(2),v(1)) - fgod(f,v(1),sin(pi*a)))
    do i = 2,N-2
      vnew(i) = v(i) - r*(fgod(f,v(j+1),v(j)) - fgod(f,v(j),v(j-1)))
    enddo
    vnew(N-1) = v(N-1) - r*(fgod(f,sin(pi*b),v(N-1)) - fgod(f,v(N-1),v(N-2)))

    v = vnew
  enddo

  do j=1,N-1
    write(unit=1, fmt=*) x(j), vnew(j)
  enddo

end program godunov
