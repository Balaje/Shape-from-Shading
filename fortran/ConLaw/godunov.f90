! Program to solve the Burger's equation using Godunov's Scheme
!     u_t + f(u)_x = 0
!       u(x,0) = { 1, x < 0
!                  0, x >= 0 }
!

program godunov

  use functions, only : f
  use flux, only : fgod

  integer, parameter :: N = 100
  real(kind = 8) :: a = -5, b = 5, h, delt, tf, t0, r, lt = 1., rt = 0.
  real(kind = 8), dimension(N-1) :: v,vnew,x
  integer :: M, i, j

  open(unit = 1, file="god.txt", status="unknown", action="write")

  h = (b-a)/N

  ! CFL Condition
  delt = 0.8*h

  ! Defining the mesh and the initial condition
  do i = 1,N-1
    x(i) = a + i*h
    if(x(i) < 0) then
      v(i) = 1.
    elseif(x(i) >= 0 .and. x(i) <= 1) then
      v(i) = 2*x(i)**3 - 3*x(i)**2 + 1
    elseif(x(i) > 1) then
      v(i) = 0.
    endif
  enddo

  tf = 1.
  t0 = 0.

  M = (tf-t0)/delt

  r = delt/h

  do j = 1,M
    vnew(1) = v(1) - r*(fgod(f,v(2),v(1)) - fgod(f,v(1),lt))
    do i = 2,N-2
      vnew(i) = v(i) - r*(fgod(f,v(j+1),v(j)) - fgod(f,v(j),v(j-1)))
    enddo
    vnew(N-1) = v(N-1) - r*(fgod(f,rt,v(N-1)) - fgod(f,v(N-1),v(N-2)))

    v = vnew
  enddo

  do j=1,N-1
    write(unit=1, fmt=*) x(j), vnew(j)
  enddo

end program godunov
