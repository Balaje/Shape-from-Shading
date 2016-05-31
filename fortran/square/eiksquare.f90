! Program to solve the eikonal equation in 2d

program eiksquare
  implicit none

  integer, parameter :: N = 200
  real (kind=8) :: a = -1., b = 1., c = -1., d = 1.
  real (kind=8) :: h, delt, dxp, dxm, dyp, dym, dx, dy, L, tol, error
  real (kind=8), dimension(N-1) :: x,y
  real (kind=8), dimension(N-1,N-1) :: u, unew
  integer :: M, i, j

  open(unit = 1, file="plot.txt", status="unknown", action="write")

  h = (b-a)/N

  delt = 0.01*h

  do j = 1,N-1
    x(j) = a + j*h
    y(j) = c + j*h
    do i=1,N-1
      u(i,j) = 0
    enddo
  enddo

  error = 100
  tol = 1e-16

  do while (error > tol)
    do j = 2,N-2
      do i = 2,N-2
        Dxp = (u(i+1,j)-u(i,j))/h
        Dxm = (u(i-1,j)-u(i,j))/h
        Dyp = (u(i,j+1)-u(i,j))/h
        Dym = (u(i,j-1)-u(i,j))/h

        Dx = min(Dxp,Dxm,0.)
        Dy = min(Dyp,Dym,0.)

        L = sqrt(Dx**2 + Dy**2) - 1.

        unew(i,j) = u(i,j) - delt*L

      enddo
    enddo
    error = maxval(abs(unew-u))
    u = unew
  enddo

  do j = 1,N-1
    do i = 1,N-1
      write(unit = 1, fmt = *) x(j), y(i), unew(i,j)
    enddo
  enddo

end program eiksquare
