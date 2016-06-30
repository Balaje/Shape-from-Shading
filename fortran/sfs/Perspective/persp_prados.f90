! Program to implement the upwinded Rouy and Tourin Scheme to
!    solve the perspective shape from shading model.

program persp1
implicit none

integer, parameter :: Nx = 201, Ny = 201
real(kind = 8), parameter :: a = 0., b = 1., c = 0., d = 1., tol = 1e-3, f = 1., pi = acos(-1.)
real(kind = 8), dimension(Nx-1,Ny-1) :: u, unew, I, exact, M, Q
real(kind = 8), dimension(Nx-1) :: x
real(kind = 8), dimension(Ny-1) :: y
real(kind = 8), dimension(3) :: Ax, Ay
real(kind = 8) :: hx, hy, delt, error, Dxp, Dxm, Dyp, Dym, Dx, Dy, L, error1, minm
integer :: j,k,Ix,Iy,p, iterations = 0

hx = (b-a)/Nx
hy = (d-c)/Ny

do j = 1,Nx-1
   x(j) = a + (j)*hx
enddo

do j = 1,Ny-1
   y(j) = c + (j)*hy
enddo

open(unit = 1, file = "blj.txt")
read(1, *) I

! Defining Parameters
minm = 0.01
do k = 1,Ny-1
  do j = 1,Nx-1
    if(I(j,k) > 0 .and. I(j,k) < minm) then
      minm = I(j,k)
    endif
  enddo
enddo

print 13, minm
13  format('Minm = ', se22.15)

do k = 1,Ny-1
   do j = 1,Nx-1
      u(j,k) = -0.5*log(minm*f**2)
      Q(j,k) = sqrt(f**2/(f**2+x(j)**2+y(k)**2))
      M(j,k) = I(j,k)*f**2/Q(j,k)
      unew(j,k) = 0
   enddo
enddo

delt = 0.5*min(hx,hy)
error = 100.

!  do p = 1,1000
do while (error> tol)
   !j = 1, k = 1
   Dxp = (u(2,1)-u(1,1))/hx
   Dxm = 0.
   Dyp = (u(1,2)-u(1,1))/hy
   Dym = 0.

   Ax = (/0.,Dxp,Dxm/)
   Ay = (/0.,Dyp,Dym/)

   Dx = max(0.,-M(2,1)*Dxp,M(1,1)*Dxm)
   Dy = max(0.,-M(1,2)*Dyp,M(1,1)*Dym)
   Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
   Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

   L = -exp(-2*u(1,1)) + sqrt(f**2*(Dx**2+Dy**2) + M(1,1)**2*(Ax(Ix)*x(1)+Ay(Iy)*y(1))**2 + M(1,1)**2*Q(1,1)**2)
   unew(1,1) = u(1,1) - delt*L

   ! j = 1, k = N+1
   Dxp = (u(2,Ny-1)-u(1,Ny-1))/hx
   Dxm = (u(1,Ny-1)-1000)/hx
   Dyp = 0.
   Dym = (u(1,Ny-1)-u(1,Ny-2))/hy

   Ax = (/0.,Dxp,Dxm/)
   Ay = (/0.,Dyp,Dym/)

   Dx = max(0.,-M(2,Ny-1)*Dxp,M(1,Ny-1)*Dxm)
   Dy = max(0.,-M(1,Ny-1)*Dyp,M(1,Ny-2)*Dym)
   Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
   Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

   L = -exp(-2*u(1,Ny-1)) + sqrt(f**2*(Dx**2+Dy**2) + M(1,Ny-1)**2*(Ax(Ix)*x(1)+Ay(Iy)*y(Ny-1))**2 + M(1,Ny-1)**2*Q(1,Ny-1)**2)
   unew(1,Ny-1) = u(1,Ny-1) - delt*L

   do k = 2,Ny-2
      do j = 2,Nx-2
        !  if(I(j,k) == 0.) then
        !     unew(j,k) = (u(j-1,k)+u(j+1,k)+u(j,k-1)+u(j,k+1))/4
        !  elseif(I(j,k) == 1.)then
        !     unew(j,k) = 0.
        !  else
            Dxp = (u(j+1,k)-u(j,k))/hx
            Dxm = (u(j,k)-u(j-1,k))/hx
            Dyp = (u(j,k+1)-u(j,k))/hy
            Dym = (u(j,k)-u(j,k-1))/hy

            Ax = (/0.,Dxp,Dxm/)
            Ay = (/0.,Dyp,Dym/)

            Dx = max(0.,-M(j+1,k)*Dxp,M(j-1,k)*Dxm)
            Dy = max(0.,-M(j,k+1)*Dyp,M(j,k-1)*Dym)
            Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
            Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

            L = -exp(-2*u(j,k)) + sqrt(f**2*(Dx**2+Dy**2) + M(j,k)**2*(Ax(Ix)*x(j)+Ay(Iy)*y(k))**2 + (M(j,k)*Q(j,k))**2 )
            unew(j,k) = u(j,k) - delt*L
        ! endif
      enddo
   enddo

   ! j = N+1, k = 1
   Dxp = 0.
   Dxm = (u(Nx-1,1)-u(Nx-2,1))/hx
   Dyp = (u(Nx-1,2)-u(Nx-1,1))/hy
   Dym = 0.

   Ax = (/0.,Dxp,Dxm/)
   Ay = (/0.,Dyp,Dym/)

   Dx = max(0.,-M(Nx-1,1)*Dxp,M(Nx-2,1)*Dxm)
   Dy = max(0.,-M(Nx-1,2)*Dyp,M(Nx-1,1)*Dym)
   Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
   Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

   L = -exp(-2*u(Nx-1,1)) + sqrt(f**2*(Dx**2+Dy**2) + M(Nx-1,1)**2*(Ax(Ix)*x(Nx-1)+Ay(Iy)*y(1))**2 + (M(Nx-1,1)*Q(Nx-1,1))**2 )
   unew(Nx-1,1) = u(Nx-1,1) - delt*L

   ! j = N+1, k = N+1
   Dxp = 0.
   Dxm = (u(Nx-1,Ny-1)-u(Nx-2,Ny-1))/hx
   Dyp = 0.
   Dym = (u(Nx-1,Ny-1)-u(Nx-1,Ny-2))/hy

   Ax = (/0.,Dxp,Dxm/)
   Ay = (/0.,Dyp,Dym/)

   Dx = max(0.,-M(Nx-1,Ny-1)*Dxp,M(Nx-2,Ny-1)*Dxm)
   Dy = max(0.,-M(Nx-1,Ny-1)*Dyp,M(Nx-1,Ny-2)*Dym)
   Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
   Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

   L = -exp(-2*u(Nx-1,Ny-1)) + sqrt(f**2*(Dx**2+Dy**2) + M(Nx-1,Ny-1)**2*(Ax(Ix)*x(Nx-1)+Ay(Iy)*y(Ny-1))**2 + (M(Nx-1,Ny-1)*Q(Nx-1,Ny-1))**2 )
   unew(Nx-1,Ny-1) = u(Nx-1,Ny-1) - delt*L

   !j = 2:N, k = 1
   do j=2,Nx-2
      Dxp = (u(j+1,1)-u(j,1))/hx
      Dxm = (u(j,1)-u(j-1,1))/hx
      Dyp = (u(j,2)-u(j,1))/hy
      Dym = 0.

      Ax = (/0.,Dxp,Dxm/)
      Ay = (/0.,Dyp,Dym/)

      Dx = max(0.,-M(j+1,1)*Dxp,M(j-1,1)*Dxm)
      Dy = max(0.,-M(j,2)*Dyp,M(j,1)*Dym)
      Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
      Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

      L = -exp(-2*u(j,1)) + sqrt(f**2*(Dx**2+Dy**2) + M(j,1)**2*(Ax(Ix)*x(j)+Ay(Iy)*y(1))**2 + (M(j,1)*Q(j,1))**2 )
      unew(j,1) = u(j,1) - delt*L
   enddo

   !j = 2:N,k=N+1
   do j=2,Nx-2
      Dxp = (u(j+1,Ny-1)-u(j,Ny-1))/hx
      Dxm = (u(j,Ny-1)-u(j-1,Ny-1))/hx
      Dyp = 0.
      Dym = (u(j,Ny-1)-u(j,Ny-2))/hy

      Ax = (/0.,Dxp,Dxm/)
      Ay = (/0.,Dyp,Dym/)

      Dx = max(0.,-M(j+1,Ny-1)*Dxp,M(j-1,Ny-1)*Dxm)
      Dy = max(0.,-M(j,Ny-1)*Dyp,M(j,Ny-2)*Dym)
      Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
      Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

      L = -exp(-2*u(j,Ny-1)) + sqrt(f**2*(Dx**2+Dy**2) + M(j,Ny-1)**2*(Ax(Ix)*x(j)+Ay(Iy)*y(Ny-1))**2 + (M(j,Ny-1)*Q(j,Ny-1))**2 )
      unew(j,Ny-1) = u(j,Ny-1) - delt*L
   enddo

   !k = 2:N, j = 1
   do k = 2,Ny-2
      Dxp = (u(2,k)-u(1,k))/hx
      Dxm = 0.
      Dyp = (u(1,k+1)-u(1,k))/hy
      Dym = (u(1,k)-u(1,k-1))/hy

      Ax = (/0.,Dxp,Dxm/)
      Ay = (/0.,Dyp,Dym/)

      Dx = max(0.,-M(2,k)*Dxp,M(1,k)*Dxm)
      Dy = max(0.,-M(1,k+1)*Dyp,M(1,k-1)*Dym)
      Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
      Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

      L = -exp(-2*u(1,k)) + sqrt(f**2*(Dx**2+Dy**2) + M(1,k)**2*(Ax(Ix)*x(1)+Ay(Iy)*y(k))**2 + (M(1,k)*Q(1,k))**2 )
      unew(1,k) = u(1,k) - delt*L
   enddo

   !k=2:N, j = N+1
   do k=2,Nx-2
      Dxp = 0.
      Dxm = (u(Nx-1,k)-u(Nx-2,k))/hx
      Dyp = (u(Nx-1,k+1)-u(Nx-1,k))/hy
      Dym = (u(Nx-1,k)-u(Nx-1,k-1))/hy

      Ax = (/0.,Dxp,Dxm/)
      Ay = (/0.,Dyp,Dym/)

      Dx = max(0.,-M(Nx-1,k)*Dxp,M(Nx-2,k)*Dxm)
      Dy = max(0.,-M(Nx-1,k+1)*Dyp,M(Nx-1,k-1)*Dym)
      Ix = maxloc((/0.,-Dxp,Dxm/),dim=1)
      Iy = maxloc((/0.,-Dyp,Dym/),dim=1)

      L = -exp(-2*u(Nx-1,k)) + sqrt(f**2*(Dx**2+Dy**2) + M(Nx-1,k)**2*(Ax(Ix)*x(Nx-1)+Ay(Iy)*y(k))**2 + (M(Nx-1,k)*Q(Nx-1,k))**2 )
      unew(Nx-1,k) = u(Nx-1,k) - delt*L
   enddo

   error = maxval(abs(u-unew))
!   print 11, error
!11 format("error=", se22.15)

 u = unew

 iterations = iterations + 1
 if(mod(iterations,100) == 0) then
   print 12, iterations
   12 format("Total Iterations = ", i10)
   print 11, error
   11 format("error=", se22.15)
 endif
enddo

open(unit = 1,file = "op_neumannf.txt")
do j=1,Nx-1
   do k=1,Ny-1
      write(1,*) x(j), y(k),-exp(unew(j,k))
   enddo
enddo



end program persp1
