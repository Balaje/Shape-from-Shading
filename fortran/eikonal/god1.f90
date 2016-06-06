  ! Fortran code to solve the 1D Eikonal equation
  ! with initial and boundary data
  !     u_t - |u_x| = -1
  !      u(x,0) = 0
  !     u(-1,t) = u(1,t) = 0

program god1
  use func, only : f1,flux1

  implicit none

  integer , parameter :: N = 100
  real(kind=8), parameter :: a = -1., b = 1., eps = 1., tol = 1e-16
  real(kind=8), dimension(N-1) :: x,u,unew
  real(kind=8) :: h,delt,tf,t0,error
  integer :: i,j,k

  h = (b-a)/N

  do i=1,N-1
     x(i) = a + i*h
     u(i) = 0.
  enddo

  delt = eps*h

  error = 100.

do while (error > tol)
     unew(1) = u(1) - delt*(flux1(f1,(u(1))/h,(u(2)-u(1))/h) + 1)
     do j = 2,N-2
        unew(j) = u(j) - delt*(flux1(f1,(u(j)-u(j-1))/h,(u(j+1)-u(j))/h) + 1)
     enddo
     unew(N-1) = u(N-1) - delt*(flux1(f1,(u(N-1)-u(N-2))/h,(-u(N-1))/h) + 1)

     error = maxval(abs(u-unew))

     print 11, error
11   format('Error = ',se22.15)

     u = unew
  enddo

  ! Write data onto a file
  open(unit = 1, file = "concave.txt", status="unknown", action="write")
  do i=1,N-1
     write(unit = 1, fmt = *)x(i),unew(i)
  enddo


end program god1
