! Program to implement the upwinded Rouy and Tourin Scheme to
!    solve the perspective shape from shading model.

program orderpersp
implicit none

integer, parameter :: q = 5
real(kind = 8), parameter :: a = 0., b = 1., c = 0., d = 1., tol = 1e-8, f = 1., pi = acos(-1.)
real(kind = 8), dimension(q) :: h, error1
real(kind = 8), dimension(q-1) :: order

real(kind = 8), dimension(3) :: Ax, Ay
real(kind = 8) :: error, Dxp, Dxm, Dyp, Dym, Dx, Dy, L, delt

real(kind = 8), dimension(:,:), allocatable :: u, unew, I, exact
real(kind = 8), dimension(:), allocatable :: x,y
integer :: N = 10,j,k,Ix,Iy,p

do p = 1,q

   h(p) = (b-a)/N;

   allocate(x(N+1))
   allocate(y(N+1))

   do j = 1,N+1
      x(j) = a + (j-1)*h(p)
      y(j) = c + (j-1)*h(p)
   enddo

!open(unit = 1, file = "mozart.txt")
!read(1, *) I

   allocate(u(N+1,N+1))
   allocate(unew(N+1,N+1))
   allocate(I(N+1,N+1))
   do k = 1,N+1
      do j = 1,N+1
         I(j,k) = f/sqrt(x(j)**2+y(k)**2+f**2)/sqrt(f**2*((2*pi*cos(2*pi*x(j))*sin(2*pi*y(k)))**2+(2*pi*sin(2*pi*x(j))*cos(2*pi*y(k)))**2) + (2*pi*x(j)*cos(2*pi*x(j))*sin(2*pi*y(k)) + 2*pi*y(k)*sin(2*pi*x(j))*cos(2*pi*y(k)))**2 + (f**2/(x(j)**2+y(k)**2+f**2)))
         u(j,k) = 0
         unew(j,k) = 0
      enddo
   enddo

   delt = h(p)/sqrt(2.)
   error = 100.

   do while (error> tol)
      do j = 2,N
         do k = 2,N
            if(abs(x(j)-0.5) < 1e-9 .and. abs(y(k)-0.5) < 1e-9) then
               unew(j,k) = 0.
            elseif(abs(x(j)-0.75) < 1e-9 .and. abs(y(k)-0.75) < 1e-9) then
               unew(j,k) = 1.
            elseif(abs(x(j)-0.25) < 1e-9 .and. abs(y(k)-0.25) < 1e-9) then
               unew(j,k) = 1.
            elseif(abs(x(j)-0.25) < 1e-9 .and. abs(y(k)-0.75) < 1e-9) then
               unew(j,k) = -1.
            elseif(abs(x(j)-0.75) < 1e-9 .and. abs(y(k)-0.25) < 1e-9) then
               unew(j,k) = -1.
            else
               Dxp = (u(j+1,k)-u(j,k))/h(p)
               Dxm = (u(j,k)-u(j-1,k))/h(p)
               Dyp = (u(j,k+1)-u(j,k))/h(p)
               Dym = (u(j,k)-u(j,k-1))/h(p)

               Ax = (/0.,Dxp,Dxm/)
               Ay = (/0.,Dyp,Dym/)

               Dx = max(0.,-Dxp,Dxm)
               Dy = max(0.,-Dyp,Dym)
               Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
               Iy = maxloc((/0.,-Dyp,Dym/),dim=1)


               L = I(j,k)*sqrt(f**2*(Dx**2+Dy**2) + (Dx*x(j)+Dy*y(k))**2 + f**2/(x(j)**2+y(k)**2+f**2)) - f/sqrt(x(j)**2+y(k)**2+f**2)

               unew(j,k) = u(j,k) - delt*L
            endif
         enddo
      enddo
      error = maxval(abs(u-unew))
!      print 11, error
!11    format("error=", se22.15)

      u = unew
   enddo

   allocate(exact(N+1,N+1))
   do j=1,N+1
      do k=1,N+1
         exact(j,k) = sin(2*pi*x(j))*sin(2*pi*y(k))
      enddo
   enddo

   error1(p) = maxval(abs(unew-exact))
   print 11, N ,error1(p)
11 format('In N =  ',i6,' error = ',se22.15)
   N = N*2


   deallocate(x)
   deallocate(y)
   deallocate(u)
   deallocate(unew)
   deallocate(I)
   deallocate(exact)

enddo

do j = 1,q-1
   order(j) = log(error1(j)/error1(j+1))/log(2.)
enddo

print 12, order
12 format('Order = ', se22.15)

end program orderpersp
