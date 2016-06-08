! Program to solve eikonal equation on a circle

program eikcircle
  implicit none
  real(kind=8),parameter::a=-1.,b=1.,c=-1.,d=1.
  integer,parameter::N = 100
  real (kind=8) :: h, delt, dxp, dxm, dyp, dym, dx, dy, L, tol, error,eps, error1, minm = 0.
  real (kind=8), dimension(N+1) :: x,y
  real (kind=8), dimension(N+1,N+1) :: u, unew, z
  real (kind=8), dimension(5) :: arr
  integer,dimension(N+1,N+1)::inout
  integer :: M, i, j, k, count = 0, t

  open(unit = 1, file="plot.txt", status="unknown", action="write")

  h = (b-a)/N

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
          minm = 0.
         if (inout(i,j) == 0) then
            minm = h
          Dxp = (u(i+1,j)-u(i,j))/h
          Dxm = (u(i-1,j)-u(i,j))/h
          Dyp = (u(i,j+1)-u(i,j))/h
          Dym = (u(i,j-1)-u(i,j))/h
          if (inout(i,j+1) == 1) then
             eps = sqrt(1- x(i)**2) - y(j);
             if(eps < minm) then
                minm = eps
             endif
             Dyp = (- u(i,j))/abs(eps);

                !BOTTOM NO
           elseif (inout(i,j-1) == 1) then 
              eps = y(j) + sqrt(1 - x(i)**2);
              if(eps < minm) then
                 minm = eps
              endif
              Dym = -(u(i,j))/abs(eps);

                !LEFT NO
           elseif (inout(i-1,j) == 1) then
              eps = x(i) + sqrt(1 - y(j)**2);
              if(eps < minm) then
                 minm = eps
              endif
              Dxm = -(u(i,j))/abs(eps);

                !RIGHT NO
           elseif (inout(i+1,j) == 1) then
              eps = - x(i) + sqrt(1 - y(j)**2);
              if(eps < minm) then
                 minm = eps
              endif
              Dxp = (- u(i,j))/abs(eps);
           endif

           delt = 0.6*minm

           Dx = Dxp
           if(Dxm < Dx) then
              Dx = Dxm
           elseif(0 < Dx) then
              Dx = 0.
           endif

           Dy = Dyp
           if(Dym < Dy) then
              Dy = Dym
           elseif(0 < Dy) then
              Dy = 0.
           endif

!          Dx = min(Dxp,Dxm,0.);
!          Dy = min(Dyp,Dym,0.);
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
11 format('Error : ',es22.15)

  do j = 1,N+1
    do i = 1,N+1
      write(unit = 1, fmt = *) x(i), y(j), unew(i,j)
    enddo
  enddo

end program eikcircle
