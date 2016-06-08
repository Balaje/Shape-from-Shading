! Program to compute the order of convergence for a 2d
! Eikonal Equation

program circleorder
  implicit none

  integer, parameter :: q = 9
  real (kind = 8), parameter :: a = -1., b = 1., c = -1., d = 1., tol = 1e-16
  real (kind = 8), dimension(q) :: h, error1
  real (kind = 8), dimension(q-1) :: order
  real (kind = 8) :: error, dxp, dxm, dym, dyp, eps, L, dx, dy, delt, minm
  integer :: N = 10, i,j,k,p

  real (kind = 8), dimension(:), allocatable :: x,y
  real (kind = 8), dimension(:,:), allocatable :: u, unew, uexact
  integer, dimension(:,:), allocatable :: inout

  do p = 1,q
    h(p) = (b-a)/N

    allocate(x(N+1))
    allocate(y(N+1))

    do j = 1,N+1
      x(j) = a + (j-1)*h(p)
      y(j) = c + (j-1)*h(p)
    enddo

    allocate(u(N+1,N+1))
    allocate(unew(N+1,N+1))
    do j = 1,N+1
      do i = 1,N+1
        u(i,j) = 0.
        unew(i,j) = 0.
      enddo
    enddo

    allocate(inout(N+1,N+1))
    do j = 1,N+1
      do i = 1,N+1
        if (x(i)**2 + y(j)**2 >= 1) then
          inout(i,j) = 1
        else
          inout(i,j) = 0
        endif
      enddo
    enddo

    error = 100

    do while (error > tol)
      do j = 2,N
        do i = 2,N
           if (inout(i,j) == 0) then
            minm = h(p)
            Dxp = (u(i+1,j)-u(i,j))/h(p)
            Dxm = (u(i-1,j)-u(i,j))/h(p)
            Dyp = (u(i,j+1)-u(i,j))/h(p)
            Dym = (u(i,j-1)-u(i,j))/h(p)
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

            L = sqrt(Dx**2+Dy**2)-1.;
            unew(i,j) = u(i,j) - delt*L;

          endif
        enddo
      enddo
      error = maxval(abs(unew-u))
      u = unew
    enddo

    allocate(uexact(N+1,N+1))
    do j = 1,N+1
      do i =1,N+1
          uexact(i,j) = 1-sqrt(x(i)**2+y(j)**2)
          if (uexact(i,j) < 0) then
              uexact(i,j) = 0
          endif
      enddo
    enddo

    error1(p) = maxval(abs(unew-uexact))
    print 11, N ,error1(p)
11    format('In N =  ',i6,' error = ',se22.15)
    N = N*2


    deallocate(x)
    deallocate(y)
    deallocate(u)
    deallocate(unew)
    deallocate(inout)
    deallocate(uexact)
  enddo

  do j = 1,q-1
    order(j) = log(error1(j)/error1(j+1))/log(2.)
  enddo

  print 12, order
12    format('Order = ', se22.15)
end program circleorder
