  ! Program to implement the upwinded Rouy and Tourin Scheme to
  !    solve the perspective shape from shading model.

program persp
  implicit none

  integer, parameter :: N = 200
  real(kind = 8), parameter :: a = 0., b = 1., c = 0., d = 1., tol = 1e-8, f = 1., pi = acos(-1.)
  real(kind = 8), dimension(N+1,N+1) :: u, unew, I, exact
  real(kind = 8), dimension(N+1) :: x,y
  real(kind = 8), dimension(3) :: Ax, Ay
  real(kind = 8) :: h, delt, error, Dxp, Dxm, Dyp, Dym, Dx, Dy, L, error1
  integer :: j,k,Ix,Iy,p

  h = (b-a)/N;

  do j = 1,N+1
     x(j) = a + (j-1)*h
     y(j) = c + (j-1)*h
  enddo

  open(unit = 1, file = "moz.txt")
  read(1, *) I

  do k = 1,N+1
     do j = 1,N+1
        !I(j,k) = f/sqrt(x(j)**2+y(k)**2+f**2)/sqrt(f**2*((2*pi*cos(2*pi*x(j))*sin(2*pi*y(k)))**2+(2*pi*sin(2*pi*x(j))*cos(2*pi*y(k)))**2) + (2*pi*x(j)*cos(2*pi*x(j))*sin(2*pi*y(k)) + 2*pi*y(k)*sin(2*pi*x(j))*cos(2*pi*y(k)))**2 + (f**2/(x(j)**2+y(k)**2+f**2)))
        !I(j,k) = 1/sqrt(2.)
        u(j,k) = -0.5*log(I(j,k)*f**2)
        unew(j,k) = 0
     enddo
  enddo

  delt = h/sqrt(2.)
  error = 100.

  !do p = 1,400
  do while (error> tol)
     do j = 2,N
        do k = 2,N
          ! if(abs(x(j)-0.5) < 1e-9 .and. abs(y(k)-0.5) < 1e-9) then
          !    unew(j,k) = 0.
          ! elseif(abs(x(j)-0.75) < 1e-9 .and. abs(y(k)-0.75) < 1e-9) then
          !    unew(j,k) = 1.
          ! elseif(abs(x(j)-0.25) < 1e-9 .and. abs(y(k)-0.25) < 1e-9) then
          !    unew(j,k) = 1.
          ! elseif(abs(x(j)-0.25) < 1e-9 .and. abs(y(k)-0.75) < 1e-9) then
          !    unew(j,k) = -1.
          ! elseif(abs(x(j)-0.75) < 1e-9 .and. abs(y(k)-0.25) < 1e-9) then
          !    unew(j,k) = -1.
           ! else
           if(I(j,k) == 0.) then
              unew(j,k) = (u(j-1,k)+u(j+1,k)+u(j,k-1)+u(j,k+1))/4
           elseif(I(j,k)==1.) then
              unew(j,k) = 0.
           else

           Dxp = (u(j+1,k)-u(j,k))/h
           Dxm = (u(j,k)-u(j-1,k))/h
           Dyp = (u(j,k+1)-u(j,k))/h
           Dym = (u(j,k)-u(j,k-1))/h

           Ax = (/0.,Dxp,Dxm/)
           Ay = (/0.,Dyp,Dym/)

           Dx = max(0.,-I(j+1,k)*Dxp,I(j-1,k)*Dxm)
           Dy = max(0.,-I(j,k+1)*Dyp,I(j,k-1)*Dym)
           Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
           Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

           L = sqrt(f**2*(Dx**2+Dy**2) + I(j,k)**2*(Ax(Ix)*x(j)+Ay(Iy)*y(k))**2 + I(j,k)**2*f**2/(x(j)**2+y(k)**2+f**2)) - f/sqrt(x(j)**2+y(k)**2+f**2)

           unew(j,k) = u(j,k) - delt*L
           endif
        enddo
     enddo
     error = maxval(abs(u-unew))
     print 11, error
11   format("error=", se22.15)

     u = unew
  enddo

  open(unit = 1,file = "op_dirich.txt")
  do j=1,N+1
     do k=1,N+1
        write(1,*) x(j), y(k), -exp(unew(j,k))
     enddo
  enddo

  error1 = maxval(abs(unew-exact))

end program persp
