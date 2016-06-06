  ! Fortran code to solve the 1-D Eikonal equation
  !     with initial data and boundary data
  !      u_t + |u_x| = 1
  !        u(x,0) = 0
  !          u(-1,t) = u(1,t) = 0

program god

  use func, only : f,flux
  implicit none

  integer, parameter :: N = 100
  real(kind = 8), parameter :: a = -1., b = 1, tol = 1e-8
  real(kind = 8) :: h,error,delt,ur=0.
  real(kind = 8), dimension(N-1) :: u, unew, x
  integer :: i,j,k,M

  h = (b-a)/N
  do i=1,N-1
     x(i) = a+(i)*h
     u(i) = 0.
  enddo

  error = 100

  delt = 1*h

  do while (error > tol)
     unew(1) = u(1) - delt*(flux(f,(u(1))/h,(u(2)-u(1))/h)-1)
     do j = 2,N-2
        unew(j) = u(j) - delt*(flux(f,(u(j)-u(j-1))/h,(u(j+1)-u(j))/h)-1)
     enddo
     unew(N-1) = u(N-1) - delt*(flux(f,(u(N-1)-u(N-2))/h,(ur-u(N-1))/h)-1)

     error = maxval(abs(u-unew))
     print 11, error
11   format('Error = ',se22.15)

     u = unew
  enddo

  ! Write data onto a file
  open(unit = 1, file = "convex.txt", status="unknown", action="write")
  do i=1,N-1
     write(unit = 1, fmt = *)x(i),unew(i)
  enddo

end program god
