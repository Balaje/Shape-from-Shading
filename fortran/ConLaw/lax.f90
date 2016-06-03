! Program to solve the Burger's equation using Lax Friedrich's Scheme
!     u_t + f(u)_x = 0
!       u(x,0) = { 1, x < 0
!                  0, x >= 0 }
!

program lax

  use functions, only : f
  use flux, only : flax

  integer, parameter :: N = 100
  real(kind = 8) :: a = -5, b = 5, h, delt, tf, t0, r, lt = 1., rt = 0.
  real(kind = 8), dimension(N-1) :: v,vnew,x
  integer :: M, i, j

  open(unit = 2, file="lax.txt", status="unknown", action="write")

  h = (b-a)/N

  ! CFL Condition
  delt = 1*h

  ! Defining the mesh and the initial condition
  do i = 1,N-1
    x(i) = a + i*h
    v(i) = sin(x(i))
  enddo

  tf = 1.
  t0 = 0.

  M = (tf-t0)/delt

  r = delt/h

  do j = 1,M
    vnew(1) = v(1) - r*(flax(f,v(1),v(2),r) - flax(f,sin(a),v(1),r))
    do i = 2,N-2
      vnew(i) = v(i) - r*(flax(f,v(i),v(i+1),r) - flax(f,v(i-1),v(i),r))
    enddo
    vnew(N-1) = v(N-1) - r*(flax(f,v(N-1),sin(b),r) - flax(f,v(N-2),v(N-1),r))

    v = vnew
  enddo

  do j=1,N-1
    write(unit=2, fmt=*) x(j), vnew(j)
  enddo

end program lax
