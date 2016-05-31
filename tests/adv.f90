! Advection equation problem in Fortran

program adv

  implicit none

  integer,parameter :: N = 100, tf = 1, t0 = 0
  real (kind=8) :: a = -5., b = 5.
  real (kind=8) :: h, delt
  integer :: M, i, j
  real (kind=8), dimension(N-1) :: u, unew, uexact, x

  open(unit = 1, file = "approx.txt", status="old",action="write")

  h = (b-a)/N

  delt = 0.9*h

  ! Initial Condition
  do i = 1,N-1
    x(i) = a+i*h
    if ( x(i) < 0 ) then
      u(i) = 1
    elseif ( x(i) >= 0 .and. x(i) <= 1 ) then
      u(i) = 2*x(i)**3 - 3*x(i)**2 + 1
    elseif ( x(i) > 1 ) then
      u(i) = 0
    end if
  enddo

  M = (tf-t0)/delt

  do i = 1,M
    unew(1) = u(1) - delt/h*(u(1) - 1)
    do j = 2,N-1
      unew(j) = u(j) - delt/h*(u(j) - u(j-1))
    enddo

    u = unew
  enddo

  do i = 1,N-1
    write (unit = 1, fmt = *) x(i),unew(i)
    if ( x(i)-tf < 0 ) then
      uexact(i) = 1
    elseif ( x(i)-tf >= 0 .and. x(i)-tf <= 1 ) then
      uexact(i) = 2*(x(i)-tf)**3 - 3*(x(i)-tf)**2 + 1
    elseif ( x(i)-tf > 1 ) then
      uexact(i) = 0
    end if
    write (unit = 1,fmt = *) x(i),uexact(i)
  enddo

end program adv
