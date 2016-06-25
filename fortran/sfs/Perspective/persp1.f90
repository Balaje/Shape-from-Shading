  ! Program to implement the upwinded Rouy and Tourin Scheme to
  !    solve the perspective shape from shading model.

program persp1
  implicit none

  integer, parameter :: N = 200
  real(kind = 8), parameter :: a = 0., b = 1., c = 0., d = 1., tol = 1e-15, f = 1., pi = acos(-1.)
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
        u(j,k) = -0.5*log(I(j,k)*f**2)
        unew(j,k) = 0
     enddo
  enddo

  delt = 0.5*h
  error = 100.

!  do p = 1,1000
  do while (error> tol)
     !j = 1, k = 1
     Dxp = (u(2,1)-u(1,1))/h
     Dxm = 0.
     Dyp = (u(1,2)-u(1,1))/h
     Dym = 0.

     Ax = (/0.,Dxp,Dxm/)
     Ay = (/0.,Dyp,Dym/)

     Dx = max(0.,-Dxp,Dxm)
     Dy = max(0.,-Dyp,Dym)
     Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
     Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

     L = I(1,1)*sqrt(f**2*(Dx**2+Dy**2) + (Ax(Ix)*x(1)+Ay(Iy)*y(1))**2 + f**2/(x(1)**2+y(1)**2+f**2)) - f/sqrt(x(1)**2+y(1)**2+f**2)
     unew(1,1) = u(1,1) - delt*L

     ! j = 1, k = N+1
     Dxp = (u(2,N+1)-u(2,N+1))/h
     Dxm = 0.
     Dyp = 0.
     Dym = (u(1,N+1)-u(1,N))/h

     Ax = (/0.,Dxp,Dxm/)
     Ay = (/0.,Dyp,Dym/)

     Dx = max(0.,-Dxp,Dxm)
     Dy = max(0.,-Dyp,Dym)
     Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
     Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

     L = I(1,N+1)*sqrt(f**2*(Dx**2+Dy**2) + (Ax(Ix)*x(1)+Ay(Iy)*y(N+1))**2 + f**2/(x(1)**2+y(N+1)**2+f**2)) - f/sqrt(x(1)**2+y(N+1)**2+f**2)
     unew(1,N+1) = u(1,N+1) - delt*L

     do k = 2,N
        do j = 2,N
           if(I(j,k) == 0) then
              unew(j,k) = (u(j-1,k)+u(j+1,k)+u(j,k+1)+u(j,k-1))/4
           elseif(I(j,k) == 1) then
              unew(j,k) = 0
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

     ! j = N+1, k = 1
     Dxp = 0.
     Dxm = (u(N+1,1)-u(N,1))/h
     Dyp = (u(N+1,2)-u(N+1,1))/h
     Dym = 0.

     Ax = (/0.,Dxp,Dxm/)
     Ay = (/0.,Dyp,Dym/)

     Dx = max(0.,-Dxp,Dxm)
     Dy = max(0.,-Dyp,Dym)
     Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
     Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

     L = I(N+1,1)*sqrt(f**2*(Dx**2+Dy**2) + (Ax(Ix)*x(N+1)+Ay(Iy)*y(1))**2 + f**2/(x(N+1)**2+y(1)**2+f**2)) - f/sqrt(x(N+1)**2+y(1)**2+f**2)
     unew(N+1,1) = u(N+1,1) - delt*L

     ! j = N+1, k = N+1
     Dxp = 0.
     Dxm = (u(N+1,N+1)-u(N,N+1))/h
     Dyp = 0.
     Dym = (u(N+1,N+1)-u(N+1,N))/h

     Ax = (/0.,Dxp,Dxm/)
     Ay = (/0.,Dyp,Dym/)

     Dx = max(0.,-Dxp,Dxm)
     Dy = max(0.,-Dyp,Dym)
     Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
     Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

     L = I(N+1,N+1)*sqrt(f**2*(Dx**2+Dy**2) + (Ax(Ix)*x(N+1)+Ay(Iy)*y(N+1))**2 + f**2/(x(N+1)**2+y(N+1)**2+f**2)) - f/sqrt(x(N+1)**2+y(N+1)**2+f**2)
     unew(N+1,N+1) = u(N+1,N+1) - delt*L

     !j = 2:N, k = 1
     do j=2,N
        Dxp = (u(j+1,1)-u(j,1))/h
        Dxm = (u(j,1)-u(j-1,1))/h
        Dyp = (u(j,2)-u(j,1))/h
        Dym = 0

        Ax = (/0.,Dxp,Dxm/)
        Ay = (/0.,Dyp,Dym/)

        Dx = max(0.,-Dxp,Dxm)
        Dy = max(0.,-Dyp,Dym)
        Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
        Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

        L = I(j,1)*sqrt(f**2*(Dx**2+Dy**2) + (Ax(Ix)*x(j)+Ay(Iy)*y(1))**2 + f**2/(x(j)**2+y(1)**2+f**2)) - f/sqrt(x(j)**2+y(1)**2+f**2)
        unew(j,1) = u(j,1) - delt*L
     enddo

     !j = 2:N,k=N+1
     do j=2,N
        Dxp = (u(j+1,N+1)-u(j,N+1))/h
        Dxm = (u(j,N+1)-u(j-1,N+1))/h
        Dyp = 0.
        Dym = (u(j,N+1)-u(j,N))/h

        Ax = (/0.,Dxp,Dxm/)
        Ay = (/0.,Dyp,Dym/)

        Dx = max(0.,-Dxp,Dxm)
        Dy = max(0.,-Dyp,Dym)
        Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
        Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

        L = I(j,N+1)*sqrt(f**2*(Dx**2+Dy**2) + (Ax(Ix)*x(j)+Ay(Iy)*y(N+1))**2 + f**2/(x(j)**2+y(N+1)**2+f**2)) - f/sqrt(x(j)**2+y(N+1)**2+f**2)
        unew(j,N+1) = u(j,N+1) - delt*L
     enddo


     error = maxval(abs(u-unew))
!     print 11, error
!11   format("error=", se22.15)

     u = unew
  enddo

  open(unit = 1,file = "op_neumann.txt")
  do j=1,N+1
     do k=1,N+1
        write(1,*) x(j), y(k),-exp(unew(j,k))
     enddo
  enddo

 end program persp1
