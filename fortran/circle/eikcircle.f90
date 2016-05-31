! Program to solve eikonal equation on a circle

program eikcircle
  implicit none
  real(kind=8),parameter::a=-1.,b=1.,c=-1.,d=1.
  integer,parameter::N = 10
  real (kind=8) :: h, delt, dxp, dxm, dyp, dym, dx, dy, L, tol, error,eps, error1
  real (kind=8), dimension(N+1) :: x,y
  real (kind=8), dimension(N+1,N+1) :: u, unew, z
  integer,dimension(N+1,N+1)::inout
  integer :: M, i, j, count = 0

  open(unit = 1, file="plot.txt", status="unknown", action="write")

  h = (b-a)/N

  delt = 0.01*h

  do i = 1,N+1
    x(i) = a + (i-1)*h
    y(i) = c + (i-1)*h
  enddo

  do i = 1,N+1
    do j=1,N+1
      u(i,j) = 0.
    enddo
  enddo

  error = 100
  tol = 1e-16

  do j = 1,N+1
    do i = 1,N+1
      if (x(i)**2 + y(j)**2 >= 1) then
        inout(i,j) = 1
      else
        inout(i,j) = 0
      endif
    enddo
  enddo

  do while (error > tol)
    do j = 2,N
      do i = 2,N
        if (inout(i,j) == 0) then
          Dxp = (u(i+1,j)-u(i,j))/h
          Dxm = (u(i-1,j)-u(i,j))/h
          Dyp = (u(i,j+1)-u(i,j))/h
          Dym = (u(i,j-1)-u(i,j))/h
          if (inout(i,j+1) == 1) then
              eps = sqrt(1- x(i)**2) - y(j);
              Dyp = (- u(i,j))/abs(eps);

                !BOTTOM NO
          elseif (inout(i,j-1) == 1) then
              eps = y(j) + sqrt(1 - x(i)**2);
              Dym = -(u(i,j))/abs(eps);

                !LEFT NO
          elseif (inout(i-1,j) == 1) then
              eps = x(i) + sqrt(1 - y(j)**2);
              Dxm = -(u(i,j))/abs(eps);

                !RIGHT NO
          elseif (inout(i+1,j) == 1) then
              eps = - x(i) + sqrt(1 - y(j)**2);
              Dxp = (- u(i,j))/abs(eps);
          endif
          Dx = min(Dxp,Dxm,0.);
          Dy = min(Dyp,Dym,0.);
          L = sqrt(Dx**2+Dy**2)-1.;
          unew(i,j) = u(i,j) - delt*L;

        endif
      enddo
    enddo
    error = maxval(abs(unew-u))
    u = unew
  enddo

  do j = 1,N+1
    do i =1,N+1
        z(i,j) = 1-sqrt(x(i)**2+y(j)**2)
        if (z(i,j) < 0) then
            z(i,j) = 0
        endif
    enddo
  enddo

  error1 = maxval(abs(unew-z))

  print 11, error1
11  format('Error : ',es22.15)

  do j = 1,N+1
    do i = 1,N+1
      write(unit = 1, fmt = *) x(i), y(j), unew(i,j)
    enddo
  enddo

end program eikcircle
